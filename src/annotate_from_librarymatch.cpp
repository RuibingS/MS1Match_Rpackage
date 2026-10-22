// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <numeric>

#include "utils_annotation.h"
#include "utils_sort.h"
#include "utils_check.h"
#include "utils_isotope.h"

using namespace Rcpp;

// ──  internal helpers ───────────────────────────────────────────────────

// trim whitespace
static inline std::string trim(const std::string& s) {
  size_t start = s.find_first_not_of(" \t\r\n");
  if (start == std::string::npos) return "";
  size_t end = s.find_last_not_of(" \t\r\n");
  return s.substr(start, end - start + 1);
}

static inline std::string tolower_copy(std::string s) {
  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  return s;
}

// OR-condition tolerance: max(abs_Da, mz * ppm * 1e-6)
static inline double mz_tol_mix(double mz, double mz_threshold, double ppm) {
  return std::max(mz_threshold, mz * ppm * 1e-6);
}

struct KeyHash {
  std::size_t operator()(const std::string& s) const noexcept {
    return std::hash<std::string>{}(s);
  }
};

// Adduct match confidence tier from ppm error
static inline std::string ppm_to_tier(double ppm_err) {
  if (ppm_err <= 5.0)  return "Level1_Confirmed";
  if (ppm_err <= 10.0) return "Level2_Probable";
  if (ppm_err <= 20.0) return "Level3_Tentative";
  return "Level4_Questionable";   // 20-30 ppm; >30 rejected before this is called
}

// Pairwise neutral-loss annotation on a mz-sorted feature set
// Helper: extract the global numeric iso group index from a label string
// Handles both in-group labels ("[grpN_k][M...]") and global labels ("[k][M...]")
static std::string extract_iso_gidx(const std::string& lbl) {
  size_t b1 = lbl.find('[');
  size_t e1 = lbl.find(']');
  if (b1 == std::string::npos || e1 == std::string::npos) return "";
  std::string inner = lbl.substr(b1 + 1, e1 - b1 - 1);
  // Strip "grpN_" prefix: "grp9_447" to "447"
  size_t us = inner.rfind('_');
  if (us != std::string::npos) {
    // Only strip if the part before '_' starts with "grp"
    std::string prefix = inner.substr(0, us);
    if (prefix.size() >= 3 &&
        prefix[0] == 'g' && prefix[1] == 'r' && prefix[2] == 'p')
      return inner.substr(us + 1);
  }
  return inner;
}

// Helper: return true if an iso label string
static bool iso_label_is_satellite(const std::string& lbl) {
  if (lbl.empty()) return false;
  // Strip _exp suffix for parsing
  std::string l = lbl;
  size_t ep = l.find("_exp");
  if (ep != std::string::npos) l = l.substr(0, ep);
  size_t e1 = l.find(']');
  if (e1 == std::string::npos) return false;
  size_t b2 = l.find('[', e1 + 1);
  if (b2 == std::string::npos) return false;
  std::string ion = l.substr(b2 + 1);  // e.g. "M+1]2+" or "M]+"
  return (ion.size() >= 2 && ion[0] == 'M' && ion[1] == '+');
}

static std::vector<std::string> run_nl(
    const Rcpp::NumericVector&   mz_s,
    const Rcpp::NumericVector&   rt_s,
    const Rcpp::CharacterVector& fn_s,
    const Rcpp::NumericVector&   nl_mass,
    const Rcpp::CharacterVector& nl_name,
    const Rcpp::CharacterVector& nl_pos,
    const Rcpp::CharacterVector& nl_neg,
    const std::string& polarity,
    double mz_threshold,
    double nl_ppm,
    double rt_threshold,
    const std::vector<std::string>& iso_labels = std::vector<std::string>())

{
  int n = mz_s.size();
  std::vector<std::string> out(n, "");
  for (int i = 0; i < n; ++i) {
    std::vector<std::string> hits;
    for (int j = i + 1; j < n; ++j) {
      if (!iso_labels.empty() && (int)iso_labels.size() > j &&
          iso_label_is_satellite(iso_labels[j]))
        continue;

      if (!Rcpp::NumericVector::is_na(rt_s[i]) &&
          !Rcpp::NumericVector::is_na(rt_s[j])) {
          if (std::fabs(rt_s[i] - rt_s[j]) > rt_threshold) continue;
      }
      double diff = mz_s[j] - mz_s[i];
      for (int k = 0; k < nl_mass.size(); ++k) {
        bool pol_ok = (polarity == "pos")
        ? (Rcpp::as<std::string>(nl_pos[k]) == "+")
          : (Rcpp::as<std::string>(nl_neg[k]) == "+");
        if (!pol_ok) continue;
        double tol = mz_tol_mix(diff, mz_threshold, nl_ppm);
        if (std::fabs(diff - nl_mass[k]) <= tol)
          hits.push_back("<- " + Rcpp::as<std::string>(nl_name[k]) +
            " (" + Rcpp::as<std::string>(fn_s[j]) + ")");
      }
    }
    out[i] = join_annotations(hits);
  }
  return out;
}

// build a polarity-filtered, deduplicated adduct table
// fills filtered_ions and filtered_formulas in-place.
static void build_adduct_table(
    const Rcpp::DataFrame& adducts_df,
    const std::string& polarity,
    std::vector<std::string>& filtered_ions,
    std::vector<std::string>& filtered_formulas)
{
  Rcpp::CharacterVector adduct_ions = adducts_df["Ion.name"];
  Rcpp::CharacterVector adduct_pol  = adducts_df["Ion.polarity"];
  Rcpp::CharacterVector adduct_mass = adducts_df["Ion.mass"];

  std::unordered_set<std::string, KeyHash> seen;
  for (int i = 0; i < adduct_ions.size(); ++i) {
    std::string pol = tolower_copy(Rcpp::as<std::string>(adduct_pol[i]));
    if      (pol == "positive" || pol == "+" || pol == "pos") pol = "pos";
    else if (pol == "negative" || pol == "-" || pol == "neg") pol = "neg";
    if (pol != polarity) continue;

    std::string ion  = Rcpp::as<std::string>(adduct_ions[i]);
    std::string mass = Rcpp::as<std::string>(adduct_mass[i]);
    if (seen.insert(ion + "||" + mass).second) {
      filtered_ions.push_back(ion);
      filtered_formulas.push_back(mass);
    }
  }
}


// ═══════════════════════════════════════════════════════════════════════════
// Function 1 — Phase 1: global isotope run + per-group annotation
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
Rcpp::List annotate_librarymatch_phase1_cpp(
    Rcpp::List      library_matches,
    Rcpp::DataFrame FT,           // full FT
    Rcpp::DataFrame adducts_df,
    Rcpp::DataFrame neutral_loss_table,
    std::string     polarity,
    double mz_threshold    = 0.005,
    double ppm             = 30,
    double pi_ppm          = 10,
    double nl_ppm          = 5,
    double iso_ppm         = 5.0,
    double rt_threshold    = 3.0,
    double pi_rt_threshold  = 30.0,
    double rt_exp_threshold = 5.0,
    double iso_threshold   = 1.003355,
    int    max_charge    = 3,
    int    max_isotope   = 4,
    double mzabs         = 0.002)
{

  // ── normalise polarity ────────────────────────────────────────────────────
  polarity = tolower_copy(polarity);
  if      (polarity == "positive" || polarity == "+" || polarity == "pos") polarity = "pos";
  else if (polarity == "negative" || polarity == "-" || polarity == "neg") polarity = "neg";

  // ── adduct table ──────────────────────────────────────────────────────────
  std::vector<std::string> filtered_ions, filtered_formulas;
  build_adduct_table(adducts_df, polarity, filtered_ions, filtered_formulas);

  // ── neutral loss table ────────────────────────────────────────────────────
  Rcpp::NumericVector   nl_mass = neutral_loss_table["Accurate.Mass"];
  Rcpp::CharacterVector nl_name = neutral_loss_table.containsElementNamed("Formula")
    ? (Rcpp::CharacterVector)neutral_loss_table["Formula"]
  : (Rcpp::CharacterVector)neutral_loss_table["Neutral.Loss"];
  Rcpp::CharacterVector nl_pos  = neutral_loss_table["Pos"];
  Rcpp::CharacterVector nl_neg  = neutral_loss_table["Neg"];

  // ── GLOBAL ISOTOPE RUN ───────────────────────────────────────────────────
  std::unordered_map<std::string, std::string, KeyHash> iso_label_map;
  // all_feat_map: feature_name to {mz, rt, intensity}
  struct FeatCoord { double mz, rt, intensity; };
  std::unordered_map<std::string, FeatCoord, KeyHash> all_feat_map;
  {
    std::unordered_map<std::string, int, KeyHash> feat_index;
    std::vector<std::string> all_fn;
    std::vector<double>      all_mz, all_rt, all_int;

    auto add_feature = [&](const std::string& nm, double mz, double rt, double inten) {
      std::string nm_t = trim(nm);
      if (feat_index.find(nm_t) == feat_index.end()) {
        feat_index[nm_t] = (int)all_fn.size();
        all_fn.push_back(nm_t);
        all_mz.push_back(mz);
        all_rt.push_back(rt);
        all_int.push_back(inten);
      }
    };

    // (a) library-match group features
    for (int g = 0; g < library_matches.size(); ++g) {
      Rcpp::List grp = library_matches[g];
      if (!grp.containsElementNamed("Matched_features_to_peaks")) continue;
      Rcpp::CharacterVector fn_g = grp["Matched_features_to_peaks"];
      Rcpp::NumericVector   mz_g = grp.containsElementNamed("Matched_features_mz")
        ? (Rcpp::NumericVector)grp["Matched_features_mz"] : Rcpp::NumericVector();
      Rcpp::NumericVector   rt_g = grp.containsElementNamed("Matched_features_rt")
        ? (Rcpp::NumericVector)grp["Matched_features_rt"]
      : Rcpp::NumericVector(mz_g.size(), NA_REAL);
      Rcpp::NumericVector   in_g = grp.containsElementNamed("Matched_features_int")
        ? (Rcpp::NumericVector)grp["Matched_features_int"]
      : Rcpp::NumericVector(mz_g.size(), NA_REAL);
      for (int k = 0; k < fn_g.size(); ++k)
        add_feature(trim(Rcpp::as<std::string>(fn_g[k])),
                    k < mz_g.size() ? (double)mz_g[k] : NA_REAL,
                        k < rt_g.size() ? (double)rt_g[k] : NA_REAL,
                            k < in_g.size() ? (double)in_g[k] : NA_REAL);
    }

    // (b) features in FT
    if (FT.nrows() > 0 && FT.containsElementNamed("feature_name")) {
      Rcpp::CharacterVector pfn = FT["feature_name"];
      auto get_col = [&](const char* a, const char* b,
                         const char* c) -> Rcpp::NumericVector {
                           if (FT.containsElementNamed(a)) return FT[a];
                           if (FT.containsElementNamed(b)) return FT[b];
                           if (FT.containsElementNamed(c)) return FT[c];
                           return Rcpp::NumericVector(pfn.size(), NA_REAL);
                         };
      Rcpp::NumericVector pmz = get_col("mzmed", "mz",    "m/z");
      Rcpp::NumericVector prt = get_col("rtmed", "rt",    "retention_time");
      Rcpp::NumericVector pin = [&]() -> Rcpp::NumericVector {
        for (const char* nm : {"maxo","into","intb","intensity","intensity_median",
             "median_intensity","maxInt","intf"}) {
          if (FT.containsElementNamed(nm)) return FT[nm];
        }
        return Rcpp::NumericVector(pfn.size(), NA_REAL);
      }();
      for (int k = 0; k < pfn.size(); ++k)
        add_feature(trim(Rcpp::as<std::string>(pfn[k])),
                    pmz[k], prt[k], pin[k]);
    }

    if (!all_fn.empty()) {
      int iso_idx = 1;
      Rcpp::CharacterVector fn_r(all_fn.begin(), all_fn.end());
      Rcpp::NumericVector   mz_r(all_mz.begin(), all_mz.end());
      Rcpp::NumericVector   rt_r(all_rt.begin(), all_rt.end());
      Rcpp::NumericVector   in_r(all_int.begin(), all_int.end());
      Rcpp::CharacterVector iso_global = Check_Iso_Charge(
        fn_r, mz_r, rt_r, in_r,
        iso_threshold, max_charge, max_isotope,
        iso_idx, iso_ppm, mzabs, rt_threshold);
      for (size_t k = 0; k < all_fn.size(); ++k) {
        iso_label_map[all_fn[k]] = Rcpp::as<std::string>(iso_global[k]);
        all_feat_map[all_fn[k]]  = {all_mz[k], all_rt[k], all_int[k]};
      }
    }
  }

  // ── ISO group  ─────────────────────────────────────────────────

  {
    struct IsoFeat {
      std::string fn, gidx, label;
      double mz;
    };
    std::vector<IsoFeat> labeled;
    labeled.reserve(iso_label_map.size());
    for (auto& kv : iso_label_map) {
      if (kv.second.empty()) continue;
      size_t b1 = kv.second.find('['), e1 = kv.second.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) continue;
      std::string gidx = kv.second.substr(b1 + 1, e1 - b1 - 1);
      auto ait = all_feat_map.find(kv.first);
      if (ait == all_feat_map.end()) continue;
      labeled.push_back({kv.first, gidx, kv.second, ait->second.mz});
    }

    std::unordered_map<std::string, std::string, KeyHash> parent;
    std::function<std::string(const std::string&)> find =
      [&](const std::string& x) -> std::string {
        if (!parent.count(x)) return x;
        return parent[x] = find(parent[x]);
      };
    auto unite = [&](const std::string& a, const std::string& b) {
      std::string ra = find(a), rb = find(b);
      if (ra == rb) return;
      // keep smller gidx as root
      int ia = 0, ib = 0;
      try { ia = std::stoi(ra); } catch(...) {}
      try { ib = std::stoi(rb); } catch(...) {}
      if (ia <= ib) parent[rb] = ra;
      else          parent[ra] = rb;
    };

    // For each pair of labeled features, check
    // k × iso_threshold / c for any k in [1, max_isotope], c in [1, max_charge], RT gate
    int nl = (int)labeled.size();
    double mz_tol = mzabs * 2.0;  // generous tolerance for chain links
    for (int i = 0; i < nl; ++i) {
      for (int j = i + 1; j < nl; ++j) {
        if (labeled[i].gidx == labeled[j].gidx) continue; // already same group
        // RT gate: skip if RT gap > rt_threshold (different elution windows)
        auto ai = all_feat_map.find(labeled[i].fn);
        auto aj = all_feat_map.find(labeled[j].fn);
        if (ai != all_feat_map.end() && aj != all_feat_map.end()) {
          double rt_i = ai->second.rt, rt_j = aj->second.rt;
          if (!std::isnan(rt_i) && !std::isnan(rt_j) &&
              std::fabs(rt_i - rt_j) > rt_threshold) continue;
        }
        double diff = std::fabs(labeled[i].mz - labeled[j].mz);
        if (diff < 1e-6) continue;
        bool linked = false;
        for (int c = 1; c <= max_charge && !linked; ++c) {
          for (int m = 1; m <= max_isotope && !linked; ++m) {
            double expected = m * iso_threshold / (double)c;
            if (std::fabs(diff - expected) <= mz_tol) linked = true;
          }
        }
        if (linked) unite(labeled[i].gidx, labeled[j].gidx);
      }
    }

    // update iso_label_map labels
    for (auto& kv : iso_label_map) {
      if (kv.second.empty()) continue;
      size_t b1 = kv.second.find('['), e1 = kv.second.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) continue;
      std::string old_gidx = kv.second.substr(b1 + 1, e1 - b1 - 1);
      std::string new_gidx = find(old_gidx);
      if (new_gidx != old_gidx) {
        kv.second = "[" + new_gidx + "]" + kv.second.substr(e1 + 1);
      }
    }
  }

  struct IsoPoolFeat {
    std::string fn, label;
    double mz, rt, intensity;
  };
  std::unordered_map<std::string, std::vector<IsoPoolFeat>, KeyHash> iso_group_map;
  {
    for (auto& kv : iso_label_map) {
      const std::string& lbl = kv.second;
      if (lbl.empty()) continue;
      size_t b1 = lbl.find('[');
      size_t e1 = lbl.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) continue;
      std::string grp_idx = lbl.substr(b1 + 1, e1 - b1 - 1);
      auto ait = all_feat_map.find(kv.first);
      double fmz    = (ait != all_feat_map.end()) ? ait->second.mz        : NA_REAL;
      double frt    = (ait != all_feat_map.end()) ? ait->second.rt        : NA_REAL;
      double finten = (ait != all_feat_map.end()) ? ait->second.intensity : NA_REAL;
      iso_group_map[grp_idx].push_back({kv.first, lbl, fmz, frt, finten});
    }
  }

  // ── Phase 1 per-group annotation loop ────────────────────────────────────
  std::vector<Rcpp::List> rows;
  rows.reserve(library_matches.size() * 6);

  // confirmed_M output for Phase 2
  std::vector<double>      cm_M;
  std::vector<std::string> cm_name, cm_lib, cm_phase1_features;
  std::vector<double>      cm_ref_rt;
  // All feature names emitted in Phase 1 — returned so Phase 2 can skip them.
  std::unordered_set<std::string, KeyHash> phase1_fn_set;
  std::vector<std::string>                 phase1_fn_vec;
  std::unordered_map<std::string, int, KeyHash> cm_found; // key = "M_rtXXX" per RT cluster


  int iso_idx_grp = 1;

  const int n_groups = library_matches.size();

  for (int g = 0; g < n_groups; ++g) {
    Rcpp::List grp = library_matches[g];

    std::string compound_name = grp.containsElementNamed("Compound_name")
      ? Rcpp::as<std::string>(grp["Compound_name"]) : "";
    std::string lib_group = grp.containsElementNamed("Lib_group")
      ? Rcpp::as<std::string>(grp["Lib_group"]) : "";
    std::string ion_type = grp.containsElementNamed("Ion_type")
      ? Rcpp::as<std::string>(grp["Ion_type"]) : "";

    double M = grp.containsElementNamed("Compound_M")
      ? Rcpp::as<double>(grp["Compound_M"]) : NA_REAL;
    double precursor_mz = grp.containsElementNamed("Precursor_mz")
      ? Rcpp::as<double>(grp["Precursor_mz"]) : NA_REAL;
    double lib_rt = NA_REAL;
    for (const char* nm : {"Lib_rt", "RetTime", "RetentionTime",
         "RT", "rt", "retention_time"}) {
      if (grp.containsElementNamed(nm)) {
        try { lib_rt = Rcpp::as<double>(grp[nm]); } catch(...) {}
        break;
      }
    }

    Rcpp::CharacterVector fn = grp.containsElementNamed("Matched_features_to_peaks")
      ? (Rcpp::CharacterVector)grp["Matched_features_to_peaks"]
    : Rcpp::CharacterVector();
    Rcpp::NumericVector mz = grp.containsElementNamed("Matched_features_mz")
      ? (Rcpp::NumericVector)grp["Matched_features_mz"]
    : Rcpp::NumericVector();
    Rcpp::NumericVector rt = grp.containsElementNamed("Matched_features_rt")
      ? (Rcpp::NumericVector)grp["Matched_features_rt"]
    : Rcpp::NumericVector(mz.size(), NA_REAL);
    Rcpp::NumericVector inten = grp.containsElementNamed("Matched_features_int")
      ? (Rcpp::NumericVector)grp["Matched_features_int"]
    : Rcpp::NumericVector(mz.size(), NA_REAL);

    if (fn.size() == 0 || mz.size() == 0) continue;
    if (mz.size() != fn.size())
      Rcpp::stop("Group %d: Matched_features_mz length != Matched_features_to_peaks.", g);
    if (rt.size()    != mz.size()) rt    = Rcpp::NumericVector(mz.size(), NA_REAL);
    if (inten.size() != mz.size()) inten = Rcpp::NumericVector(mz.size(), NA_REAL);

    std::string pi_fn = grp.containsElementNamed("feature_name")
      ? Rcpp::as<std::string>(grp["feature_name"]) : "";
    double pi_mz = grp.containsElementNamed("feature_mz")
      ? Rcpp::as<double>(grp["feature_mz"]) : NA_REAL;
    double pi_rt = grp.containsElementNamed("feature_rt")
      ? Rcpp::as<double>(grp["feature_rt"]) : NA_REAL;
    double pi_int = grp.containsElementNamed("feature_int")
      ? Rcpp::as<double>(grp["feature_int"]) : NA_REAL;

    if (!pi_fn.empty() && !Rcpp::NumericVector::is_na(pi_mz)) {
      // check if already present
      bool already_in_fn = false;
      for (int i = 0; i < (int)fn.size(); ++i) {
        if (Rcpp::as<std::string>(fn[i]) == pi_fn) { already_in_fn = true; break; }
      }
      if (!already_in_fn) {
        // prepend PI feature to fn / mz / rt / inten
        int old_n = fn.size();
        Rcpp::CharacterVector fn2(old_n + 1);
        Rcpp::NumericVector   mz2(old_n + 1), rt2(old_n + 1), in2(old_n + 1);
        fn2[0] = pi_fn;  mz2[0] = pi_mz;
        rt2[0] = pi_rt;  in2[0] = Rcpp::NumericVector::is_na(pi_int) ? NA_REAL : pi_int;
        for (int i = 0; i < old_n; ++i) {
          fn2[i+1] = fn[i]; mz2[i+1] = mz[i];
          rt2[i+1] = rt[i]; in2[i+1] = inten[i];
        }
        fn = fn2; mz = mz2; rt = rt2; inten = in2;
      }
    }

    double inner_cos = grp.containsElementNamed("inner_cosine")
      ? Rcpp::as<double>(grp["inner_cosine"]) : NA_REAL;
    double outer_cos = grp.containsElementNamed("outer_cosine")
      ? Rcpp::as<double>(grp["outer_cosine"]) : NA_REAL;
    double match_sc  = grp.containsElementNamed("Match_peak_score")
      ? Rcpp::as<double>(grp["Match_peak_score"]) : NA_REAL;

    // theoretical adduct m/z for this M
    std::vector<double> theo_adduct_mz;
    if (!Rcpp::NumericVector::is_na(M) && !filtered_formulas.empty())
      theo_adduct_mz = parse_adduct_vector(M, filtered_formulas);

    // PI/precursor target
    double PI_target = NA_REAL;
    std::string PI_tag = "";
    {
      std::string it = tolower_copy(ion_type);
      if (!Rcpp::NumericVector::is_na(precursor_mz) &&
          it.find("precursor") != std::string::npos) {
        PI_target = precursor_mz;
        PI_tag    = (polarity == "pos") ? " [M+H]+ " : " [M-H]- ";
      } else if (!Rcpp::NumericVector::is_na(M)) {
        if (it.find("m+h") != std::string::npos) {
          PI_target = M + 1.007276; PI_tag = " [M+H]+ ";
        } else if (it.find("m-h") != std::string::npos) {
          PI_target = M - 1.007276; PI_tag = " [M-H]- ";
        }
      }
    }
    int pi_idx = -1;
    if (!Rcpp::NumericVector::is_na(PI_target)) {

      double pi_tol = mz_tol_mix(PI_target, mz_threshold, pi_ppm);
      double best = 1e99;
      for (int i = 0; i < mz.size(); ++i) {
        if (!Rcpp::NumericVector::is_na(lib_rt) &&
            !Rcpp::NumericVector::is_na(rt[i])) {
            if (std::fabs(rt[i] - lib_rt) > pi_rt_threshold) continue;
        }
        double d = std::fabs(mz[i] - PI_target);
        if (d < best) { best = d; pi_idx = i; }
      }
      if (best > pi_tol) pi_idx = -1;
    }

    // ── Sort group features by mz for NL pass and per-group iso run ──────
    int n_feat = mz.size();
    std::vector<int> ord(n_feat);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](int a, int b){ return mz[a] < mz[b]; });

    Rcpp::CharacterVector fn_s(n_feat);
    Rcpp::NumericVector   mz_s(n_feat), rt_s(n_feat), in_s(n_feat);
    for (int k = 0; k < n_feat; ++k) {
      fn_s[k] = fn[ord[k]]; mz_s[k] = mz[ord[k]];
      rt_s[k] = rt[ord[k]]; in_s[k] = inten[ord[k]];
    }

    // ── Isotope annotation: local group run (Phase 1) + global fallback ──────

    std::unordered_set<std::string, KeyHash> grp_fn_set;
    for (int k = 0; k < n_feat; ++k)
      grp_fn_set.insert(trim(Rcpp::as<std::string>(fn_s[k])));

    std::vector<std::string> iso_s(n_feat, "");
    if (n_feat >= 2) {

      Rcpp::CharacterVector loc_iso = Check_Iso_Charge(
        fn_s, mz_s, rt_s, in_s,
        iso_threshold, max_charge, max_isotope,
        iso_idx_grp, iso_ppm, mzabs, rt_threshold);

      for (int k = 0; k < n_feat; ++k) {
        std::string lbl = Rcpp::as<std::string>(loc_iso[k]);
        if (!lbl.empty()) {
          size_t b1 = lbl.find('['), e1 = lbl.find(']');
          if (b1 != std::string::npos && e1 != std::string::npos) {
            std::string inner = lbl.substr(b1 + 1, e1 - b1 - 1);
            lbl = "[grp" + std::to_string(g + 1) + "_" + inner + "]" + lbl.substr(e1 + 1);
          }
          iso_s[k] = lbl;  // Tier1b_Isotope — both partners confirmed in this group
        }
      }
    }


    if (!Rcpp::NumericVector::is_na(M) && !theo_adduct_mz.empty()) {
      size_t L = std::min(theo_adduct_mz.size(), filtered_ions.size());
      double mz_tol_strip = mz_threshold + M * 30.0 * 1e-6; // 30 ppm loose gate

      // Build map: group_index_string -> monoisotope mz
      // iso_s labels have the form "[grp<G>_<idx>][M]+" (monoisotope)
      // or "[grp<G>_<idx>][M+k]+" (satellite)
      std::unordered_map<std::string, double> grp_mono_mz;
      for (int k = 0; k < n_feat; ++k) {
        if (iso_s[k].empty()) continue;
        // Check for monoisotope pattern: ends with "][M]+" or "][M]2+" etc.
        // i.e. contains "][M]" (not "][M+")
        size_t mono_pos = iso_s[k].find("][M]");
        if (mono_pos == std::string::npos) continue;
        // Extract group key: everything up to and including the first ']'
        size_t e1 = iso_s[k].find(']');
        if (e1 == std::string::npos) continue;
        std::string grp_key = iso_s[k].substr(0, e1 + 1);
        grp_mono_mz[grp_key] = mz_s[k];
      }

      // Helper lambda: find multimer order from best adduct match for a given mz
      auto get_multimer_order = [&](double query_mz) -> double {
        double best_diff = 1e9;
        int    best_idx  = -1;
        for (size_t ai = 0; ai < L; ++ai) {
          double theo = theo_adduct_mz[ai];
          if (std::isnan(theo)) continue;
          double diff = std::fabs(query_mz - theo);
          if (diff < best_diff) { best_diff = diff; best_idx = (int)ai; }
        }
        if (best_idx < 0 || best_diff > mz_tol_strip) return 1.0;
        const std::string& ion = filtered_ions[(size_t)best_idx];
        size_t bpos = ion.find('[');
        if (bpos == std::string::npos || bpos + 1 >= ion.size() ||
            !std::isdigit(static_cast<unsigned char>(ion[bpos + 1]))) return 1.0;
        size_t dpos = bpos + 1;
        std::string digits;
        while (dpos < ion.size() &&
               std::isdigit(static_cast<unsigned char>(ion[dpos])))
          digits += ion[dpos++];
        return digits.empty() ? 1.0 : std::stod(digits);
      };

      for (int k = 0; k < n_feat; ++k) {
        if (iso_s[k].empty()) continue;

        double lookup_mz = mz_s[k];  // default: use own mz

        // If this is a satellite (contains "][M+" but not "][M]"),
        // use the monoisotope mz of the same group for adduct lookup.
        bool is_satellite = (iso_s[k].find("][M+") != std::string::npos) &&
                            (iso_s[k].find("][M]") == std::string::npos);
        if (is_satellite) {
          size_t e1 = iso_s[k].find(']');
          if (e1 != std::string::npos) {
            std::string grp_key = iso_s[k].substr(0, e1 + 1);
            auto it = grp_mono_mz.find(grp_key);
            if (it != grp_mono_mz.end())
              lookup_mz = it->second;  // use monoisotope mz for order lookup
          }
        }

        double multimer_order = get_multimer_order(lookup_mz);
        double ref_M = M * multimer_order;
        if (mz_s[k] > ref_M * 1.5)
          iso_s[k] = "";
      }
    } else if (!Rcpp::NumericVector::is_na(M)) {
      // Fallback: no adduct table available — use original monomer check
      for (int k = 0; k < n_feat; ++k) {
        if (iso_s[k].empty()) continue;
        if (mz_s[k] > M * 1.5)
          iso_s[k] = "";
      }
    }

    // ── NL within group (mz-sorted inputs required) ───────────────────────
    // Pass iso_s so run_nl can skip confirmed satellites as NL targets.
    std::vector<std::string> nl_s = run_nl(
      mz_s, rt_s, fn_s, nl_mass, nl_name, nl_pos, nl_neg,
      polarity, mz_threshold, nl_ppm, rt_threshold, iso_s);

    // ── Adduct annotation with confidence tier ────────────────────────────
    std::vector<std::string> adduct_s(n_feat,""), adduct_ppm_s(n_feat,""),
    adduct_mda_s(n_feat,""), adduct_conf_s(n_feat,"");
    if (!Rcpp::NumericVector::is_na(M) && !theo_adduct_mz.empty()) {
      size_t L = std::min(theo_adduct_mz.size(), filtered_ions.size());
      for (int i = 0; i < n_feat; ++i) {
        std::unordered_set<std::string, KeyHash> seen;
        std::vector<std::string> h_ann, h_ppm, h_mda, h_conf;
        for (size_t k = 0; k < L; ++k) {
          double theo = theo_adduct_mz[k];
          if (std::isnan(theo)) continue;
          double tol      = mz_tol_mix(theo, mz_threshold, ppm);
          double abs_diff = std::fabs(mz_s[i] - theo);
          if (abs_diff > tol) continue;
          double ppm_err = (theo > 0.0) ? (abs_diff / theo * 1e6) : NA_REAL;
          if (!std::isnan(ppm_err) && ppm_err > 30.0) continue;

          std::ostringstream ann_oss;
          ann_oss << filtered_ions[k] << " " << std::fixed << std::setprecision(6) << M;
          std::string ann = ann_oss.str();
          if (!seen.insert(ann).second) continue;

          std::ostringstream ppm_oss, mda_oss;
          if (!std::isnan(ppm_err))
            ppm_oss << std::fixed << std::setprecision(2) << ppm_err << " ppm";
          else ppm_oss << "NA";
          mda_oss << std::fixed << std::setprecision(4) << (abs_diff * 1e3) << " mDa";

          h_ann.push_back(ann);
          h_ppm.push_back(ppm_oss.str());
          h_mda.push_back(mda_oss.str());
          h_conf.push_back(std::isnan(ppm_err) ? "Level3_Tentative" : ppm_to_tier(ppm_err));
        }
        adduct_s[i]     = join_annotations(h_ann);
        adduct_ppm_s[i] = join_annotations(h_ppm);
        adduct_mda_s[i] = join_annotations(h_mda);
        adduct_conf_s[i]= join_annotations(h_conf);
      }
    }

    // unsort — map sorted positions back to original feature order
    std::vector<int> rev_ord(n_feat);
    for (int k = 0; k < n_feat; ++k) rev_ord[ord[k]] = k;

    // ── RT-cluster-based confirmed_M registration ─────────────────────────
    // Problem: when the library compound has no RT, all adduct-matched features
    // (potentially spanning hundreds of seconds) land in one flat group and
    // produce a single confirmed_M with ref_rt = NA.  Phase 2 then scans the
    // entire pool with no RT gate — producing false cross-RT matches.
    //
    // Fix: cluster the matched features by RT using single-linkage with gap
    // threshold = rt_exp_threshold (the same window used by Phase 2).
    // Each RT cluster registers its own confirmed_M entry with the median RT
    // of that cluster as ref_rt.  Phase 2 then expands independently per
    // cluster, each correctly RT-scoped.
    //
    // If the library compound DOES have an RT (lib_rt is not NA), we still
    // cluster but use lib_rt as the anchor — features outside pi_rt_threshold
    // were already excluded from the group, so clustering is a no-op in
    // the typical case.
    //
    // Features with RT = NA are collected into a separate NA-RT entry that
    // gets ref_rt = NA (Phase 2 will skip the RT gate for those, which is
    // the correct behaviour when RT is genuinely unknown).

    if (!Rcpp::NumericVector::is_na(M)) {
      // Separate features into RT-known and RT-unknown
      std::vector<std::pair<double,int>> rt_known;  // {rt, feature_index}
      std::vector<int>                   rt_na_idx;
      for (int i = 0; i < n_feat; ++i) {
        if (!Rcpp::NumericVector::is_na(rt[i]))
          rt_known.push_back({rt[i], i});
        else
          rt_na_idx.push_back(i);
      }

      // Sort RT-known features by RT for single-linkage clustering
      std::sort(rt_known.begin(), rt_known.end());

      // Single-linkage: new cluster when consecutive RT gap > rt_exp_threshold
      std::vector<std::vector<int>> clusters;
      for (auto& p : rt_known) {
        if (clusters.empty() ||
            p.first - rt_known[clusters.back().back()].first > rt_exp_threshold) {
          clusters.push_back({});
        }
        clusters.back().push_back((int)(
            &p - &rt_known[0]));  // index into rt_known
      }

      // Register one confirmed_M per RT cluster
      for (auto& cl : clusters) {
        // Compute median RT of this cluster
        std::vector<double> cl_rts;
        for (int ci : cl) cl_rts.push_back(rt_known[ci].first);
        std::sort(cl_rts.begin(), cl_rts.end());
        int mid = (int)cl_rts.size() / 2;
        double cluster_rt = (cl_rts.size() % 2 == 0)
          ? (cl_rts[mid-1] + cl_rts[mid]) / 2.0
        : cl_rts[mid];

        // Key = M + cluster_rt to allow same M at different RTs
        std::ostringstream mk;
        mk << std::fixed << std::setprecision(6) << M
           << "_rt" << std::setprecision(2) << cluster_rt;
        std::string mkey = mk.str();
        if (cm_found.find(mkey) == cm_found.end()) {
          cm_found[mkey] = (int)cm_M.size();
          cm_M.push_back(M);
          cm_name.push_back(compound_name);
          cm_lib.push_back(lib_group);
          cm_ref_rt.push_back(cluster_rt);
        }
      }

      // Register one confirmed_M for RT-unknown features (ref_rt = NA)
      // only if there are no RT-known clusters at all (otherwise Phase 2
      // will cover them via the RT-gated clusters above).
      if (clusters.empty() && !rt_na_idx.empty()) {
        std::ostringstream mk;
        mk << std::fixed << std::setprecision(6) << M << "_rtNA";
        std::string mkey = mk.str();
        if (cm_found.find(mkey) == cm_found.end()) {
          cm_found[mkey] = (int)cm_M.size();
          cm_M.push_back(M);
          cm_name.push_back(compound_name);
          cm_lib.push_back(lib_group);
          cm_ref_rt.push_back(NA_REAL);
        }
      }
    }

    // ── emit rows ─────────────────────────────────────────────────────────
    for (int i = 0; i < n_feat; ++i) {
      int si = rev_ord[i];
      std::string PI_str = "";
      if (i == pi_idx && !Rcpp::NumericVector::is_na(M)) {
        std::ostringstream ms; ms << std::fixed << std::setprecision(5) << M;
        PI_str = compound_name + PI_tag + ms.str();
      }

      // ── Confidence tier assignment ──────────────────────────────────────
      // Priority (highest → lowest):
      //   Tier1a_AnnoMatch    : adduct match within group
      //   Tier1b_Isotope        : isotope pair within group (in-group iso run)
      //   Tier1d_NeutralLoss    : neutral loss within group
      //   Tier1b_exp_Isotope    : isotope partner found in pool only (_exp)
      //   Tier2_Unexplained     : no evidence of any kind
      //
      // Note: Tier1b_exp_Isotope ranks below Tier1d_NeutralLoss because in-group
      // NL is directly linked to compound identity, whereas pool-level isotope
      // only confirms signal authenticity without identifying the compound.
      // Note: Tier1d (MS2 fragment confirmed) is not yet implemented.
      std::string tier;
      const std::string& iso_i = iso_s[si];
      bool iso_ingroup = !iso_i.empty() && iso_i.find("_exp") == std::string::npos;
      bool iso_pool    = !iso_i.empty() && iso_i.find("_exp") != std::string::npos;

      bool iso_is_anchor = false;
      if (iso_ingroup) {
        // Label format: "[grpN_k][M+j]c+" or "[grpN_k][M]c+"
        // The SECOND bracket pair is the ion type. An anchor has [M] (no +j).
        // A satellite has [M+1], [M+2], [M+1]2+, etc.
        // Find the start of the second '[' (after the first ']')
        size_t e1 = iso_i.find(']');           // end of group-index bracket
        size_t b2 = (e1 != std::string::npos)
          ? iso_i.find('[', e1 + 1)  // start of ion-type bracket
            : std::string::npos;
        if (b2 != std::string::npos) {
          // suffix of ion-type bracket, e.g. "M]+" or "M+1]+" or "M+1]2+"
          std::string ion_part = iso_i.substr(b2 + 1);  // e.g. "M]+" or "M+1]2+"
          // Anchor: starts with "M]" — no charge offset
          iso_is_anchor = (ion_part.size() >= 2 &&
            ion_part[0] == 'M' &&
            ion_part[1] == ']');
        }
      }

      // ── Suppress NL for confirmed isotope satellites ────────────────────────
      // A C13 isotope satellite ([M+1]+, [M+2]+, [M+1]2+, etc.) cannot
      // simultaneously be a genuine neutral-loss fragment of a different ion.
      // Any NL match on a satellite is a coincidental mz spacing, not a real
      // chemical relationship.  Clear nl_s[si] for all non-anchor isotope
      // features — both in-group and pool-level (_exp) labels — so that
      // spurious NL assignments are suppressed before tier resolution.
      // Monoisotope anchors ([M]+, [M]2+, [M]3+) are left untouched: a
      // monoisotopic ion can legitimately lose a neutral fragment.
      {
        auto is_satellite = [&](const std::string& lbl) -> bool {
          if (lbl.empty()) return false;
          // Strip _exp suffix for label parsing
          std::string l = lbl;
          size_t exp_pos = l.find("_exp");
          if (exp_pos != std::string::npos) l = l.substr(0, exp_pos);
          // Find second bracket pair (ion type): "[grpN_k][M+j]c+"
          size_t e1 = l.find(']');
          if (e1 == std::string::npos) return false;
          size_t b2 = l.find('[', e1 + 1);
          if (b2 == std::string::npos) return false;
          std::string ion_part = l.substr(b2 + 1); // e.g. "M+1]2+" or "M]+"
          // Satellite: starts with "M+" (has isotope offset)
          return (ion_part.size() >= 2 &&
                  ion_part[0] == 'M' &&
                  ion_part[1] == '+');
        };

        if (is_satellite(iso_s[si]))
          nl_s[si] = "";
      }

      if (!adduct_s[si].empty())             tier = "Tier1a_AnnoMatch";
      else if (iso_ingroup && !iso_is_anchor) tier = "Tier1b_Isotope";
      else if (!nl_s[si].empty())             tier = "Tier1d_NeutralLoss";
      else if (iso_pool)                      tier = "Tier1b_exp_Isotope";
      else                                    tier = "Tier2_Unexplained";

      // ── Global iso fallback for monoisotope anchors ───────────────────────
      // A Phase 1 feature may be the monoisotope of a global iso pair
      // (its satellite exists in the pool but NOT in this library group),
      // yet have an empty iso_s[si] because the local group run found no
      // partner inside the group.  Write the global iso label with "_exp"
      // suffix so the anchor's iso_annotation is populated and consistent
      // with its satellite's "[grpN_k][M+1]+_exp" label emitted by
      // Phase2_IsoExpansion.  Tier is unchanged — adduct/NL evidence takes
      // priority and the anchor is already correctly classified.
      if (iso_s[si].empty()) {
        std::string fn_str_iso = trim(Rcpp::as<std::string>(fn[i]));
        auto git = iso_label_map.find(fn_str_iso);
        if (git != iso_label_map.end() && !git->second.empty()) {
          const std::string& glbl = git->second;
          // Only write for monoisotope ([M]+ role, not a satellite [M+k]+)
          bool is_global_mono = (glbl.find("[M]") != std::string::npos &&
                                 glbl.find("[M+") == std::string::npos);
          if (is_global_mono) {
            size_t b1 = glbl.find('['), e1 = glbl.find(']');
            if (b1 != std::string::npos && e1 != std::string::npos) {
              std::string grp_num = lib_group;
              size_t us = grp_num.rfind('_');
              if (us != std::string::npos) grp_num = grp_num.substr(us + 1);

              // Prefer the index of an already-emitted locally-labeled partner
              // in the same lib_group whose mz is at iso spacing from this
              // feature (e.g. FT1476 has local label [grp7_74][M]+ from the
              // per-group run; FT1483 is its c=2 satellite and should inherit
              // index 74, not global index 194).
              std::string use_idx = glbl.substr(b1 + 1, e1 - b1 - 1); // global default
              double cur_mz_val = mz[i];
              if (!std::isnan(cur_mz_val)) {
                for (const Rcpp::List& er : rows) {
                  if (Rcpp::as<std::string>(er["Lib_group"]) != lib_group) continue;
                  std::string er_iso = Rcpp::as<std::string>(er["iso_annotation"]);
                  if (er_iso.empty()) continue;
                  // Extract index from existing label (handles both local
                  // "[grp7_74][M]+" and global-fallback "[grp7_194][M]+_exp")
                  size_t eb1 = er_iso.find('['), ee1 = er_iso.find(']');
                  if (eb1 == std::string::npos || ee1 == std::string::npos) continue;
                  std::string er_inner = er_iso.substr(eb1 + 1, ee1 - eb1 - 1);
                  // Strip grpN_ prefix
                  size_t er_us = er_inner.rfind('_');
                  std::string er_idx = (er_us != std::string::npos)
                    ? er_inner.substr(er_us + 1) : er_inner;
                  double er_mz = Rcpp::as<double>(er["mz"]);
                  if (std::isnan(er_mz)) continue;
                  double diff = std::fabs(cur_mz_val - er_mz);
                  if (diff < 1e-6) continue;
                  bool linked = false;
                  for (int c = 1; c <= max_charge && !linked; ++c) {
                    for (int m = 1; m <= max_isotope && !linked; ++m) {
                      double expected = m * iso_threshold / (double)c;
                      if (std::fabs(diff - expected) <= mzabs * 2.0) linked = true;
                    }
                  }
                  if (linked) { use_idx = er_idx; break; }
                }
              }

              iso_s[si] = "[grp" + grp_num + "_" + use_idx + "]" +
                           glbl.substr(e1 + 1) + "_exp";
            }
          }
        }
      }

      // Track every Phase 1 feature name for Phase 2 skip-set.
      {
        std::string fn_str = trim(Rcpp::as<std::string>(fn[i]));
        if (phase1_fn_set.insert(fn_str).second)
          phase1_fn_vec.push_back(fn_str);
      }
      rows.push_back(Rcpp::List::create(
          Rcpp::Named("feature_name")          = Rcpp::as<std::string>(fn[i]),
          Rcpp::Named("mz")                    = mz[i],
          Rcpp::Named("rt")                    = rt[i],
          Rcpp::Named("intensity")             = inten[i],
          Rcpp::Named("Compound_name")         = compound_name,
          Rcpp::Named("Lib_group")             = lib_group,
          Rcpp::Named("Compound_M")            = M,
          Rcpp::Named("PI_match")              = PI_str,
          Rcpp::Named("grp_anno_match")      = adduct_s[si],
          Rcpp::Named("grp_anno_ppm")        = adduct_ppm_s[si],
          Rcpp::Named("grp_anno_mda")        = adduct_mda_s[si],
          Rcpp::Named("grp_anno_confidence") = adduct_conf_s[si],
          Rcpp::Named("exp_anno_match")      = std::string(""),
          Rcpp::Named("exp_anno_ppm")        = std::string(""),
          Rcpp::Named("exp_anno_mda")        = std::string(""),
          Rcpp::Named("exp_anno_confidence") = std::string(""),
          Rcpp::Named("NL_match")              = nl_s[si],
          Rcpp::Named("iso_annotation")        = iso_s[si],
          Rcpp::Named("reference_tier")       = tier,
          Rcpp::Named("annotation_phase")      = "Phase1",
          Rcpp::Named("Match_peak_score")      = match_sc,
          Rcpp::Named("inner_cosine")          = inner_cos,
          Rcpp::Named("outer_cosine")          = outer_cos
      ));
    }
  } // end per-group loop

  // ── Phase2_IsoExpansion pass ─────────────────────────────────────────────
  // For every Phase 1 feature that received any _exp iso label (regardless of
  // whether it is [M]+, [M+1]+, [M+2]+, etc.), look up ALL other members of
  // its iso group and emit each missing partner as a standalone
  // Phase2_IsoExpansion row linked to the same compound / lib_group.
  //
  // Design decisions:
  //
  //   Trigger from ANY _exp member (not just [M]+):
  //     The old code only triggered from [M]+ rows.  This silently dropped
  //     all partner rows when the library group contained only a higher-order
  //     satellite (e.g. [M+2]+) but not the monoisotope.  Now any _exp member
  //     triggers expansion of the whole iso group.
  //
  //   Deduplicate by {gidx + lib_group}:
  //     Multiple Phase 1 features from the same iso group and same lib_group
  //     (e.g. [M+1]+_exp AND [M+2]+_exp both in Lib_group_1) must only fire
  //     the expansion once.  We track triggered {gidx||lib_group} pairs.
  //
  //   Path A — iso_group_map lookup:
  //     Emits partners that received a label from Check_Iso_Charge (i.e. they
  //     passed the mz / RT / intensity checks).  These are the "clean" cases.
  //
  //   Path B — direct pool scan:
  //     Catches partners that Check_Iso_Charge could NOT label because the
  //     triggering feature is an adduct form (e.g. [M+Na]+) whose intensity
  //     differs from the true monoisotope, causing the m=1 intensity window
  //     check to fail.  We scan the pool directly within rt_exp_threshold and
  //     a 30 mDa mz window around each expected isotope position.
  //
  //   Skip self:
  //     The triggering feature itself (pr.fn) is already in Phase 1 and must
  //     not be re-emitted.
  //
  //   phase1_fn_set updated:
  //     All emitted partners are added so Phase 2 M-expansion skips them.
  {
    // ── Pool arrays for Path B ─────────────────────────────────────────────
    std::vector<std::string> pool_fn_vec;
    std::vector<double>      pool_mz_vec, pool_rt_vec, pool_int_vec;
    if (FT.nrows() > 0 && FT.containsElementNamed("feature_name")) {
      Rcpp::CharacterVector _pfn = FT["feature_name"];
      auto _get = [&](const char* a, const char* b,
                      const char* c) -> Rcpp::NumericVector {
                        if (FT.containsElementNamed(a)) return FT[a];
                        if (FT.containsElementNamed(b)) return FT[b];
                        if (FT.containsElementNamed(c)) return FT[c];
                        return Rcpp::NumericVector(_pfn.size(), NA_REAL);
                      };
      Rcpp::NumericVector _pmz = _get("mzmed", "mz",  "m/z");
      Rcpp::NumericVector _prt = _get("rtmed", "rt",  "retention_time");
      Rcpp::NumericVector _pin = [&]() -> Rcpp::NumericVector {
        for (const char* nm : {"maxo","into","intb","intensity",
             "intensity_median","median_intensity","maxInt","intf"})
          if (FT.containsElementNamed(nm)) return FT[nm];
        return Rcpp::NumericVector(_pfn.size(), NA_REAL);
      }();
      int _np = _pfn.size();
      pool_fn_vec.reserve(_np); pool_mz_vec.reserve(_np);
      pool_rt_vec.reserve(_np); pool_int_vec.reserve(_np);
      for (int i = 0; i < _np; ++i) {
        pool_fn_vec.push_back(trim(Rcpp::as<std::string>(_pfn[i])));
        pool_mz_vec.push_back(_pmz[i]);
        pool_rt_vec.push_back(_prt[i]);
        pool_int_vec.push_back(_pin[i]);
      }
    }
    int n_pool = (int)pool_fn_vec.size();

    // ── Collect triggering Phase 1 rows ───────────────────────────────────
    // Any Phase 1 row with an _exp iso label is a potential trigger.
    // Store the best representative per {gidx + lib_group} pair so we always
    // use a consistent compound context for the emitted rows.
    struct P1IsoRow {
      std::string fn, compound_name, lib_group, iso_label, gidx;
      std::string local_idx;  // index to use in emitted labels: local in-group index
                              // when the trigger has a local Phase 1 label (e.g. "66"
                              // from "[grp7_66][M]+"), global gidx otherwise (e.g. "10").
                              // Ensures emitted satellite labels share the same index as
                              // the anchor's existing Phase 1 iso_annotation.
      double M, match_sc, inner_cos, outer_cos, mz, rt;
      double ref_rt;  // RT of the [M+H]+ anchor for this lib_group (from confirmed_M);
                      // used for pool RT gate instead of the satellite's own RT.
    };

    // Helper: extract the numeric index from a local in-group iso label.
    // "[grpN_IDX][M]+" → "IDX".  Returns "" if not a local label.
    auto extract_local_idx = [](const std::string& lbl) -> std::string {
      if (lbl.empty() || lbl.find("_exp") != std::string::npos) return "";
      size_t b1 = lbl.find('['), e1 = lbl.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) return "";
      std::string inner = lbl.substr(b1 + 1, e1 - b1 - 1);
      // Must have grpN_ prefix
      size_t us = inner.rfind('_');
      if (us == std::string::npos) return "";
      std::string prefix = inner.substr(0, us);
      if (prefix.size() < 3 || prefix[0]!='g' || prefix[1]!='r' || prefix[2]!='p') return "";
      return inner.substr(us + 1);
    };

    // Map {lib_group} → ref_rt for [M+H]+ anchor.
    // Used as the RT gate for pool scan (Bug 2 fix): isotope expansion should
    // check pool features against the [M+H]+ RT, not the satellite's own RT.
    // When a lib_group has multiple RT clusters (same compound, multiple peaks),
    // we store the one whose ref_rt is closest to the triggering satellite — this
    // is resolved per-trigger below.  For now build a simple first-hit map; the
    // per-trigger refinement happens when constructing P1IsoRow.
    std::unordered_map<std::string, double, KeyHash> libgrp_refrt_map;
    for (int ci = 0; ci < (int)cm_lib.size(); ++ci) {
      // Only insert once; for multi-cluster compounds the first cluster's RT
      // is a reasonable fallback — exact match is done per-trigger below.
      if (libgrp_refrt_map.find(cm_lib[ci]) == libgrp_refrt_map.end())
        libgrp_refrt_map[cm_lib[ci]] = cm_ref_rt[ci];
    }

    // Map {gidx||lib_group} → best trigger row (prefer [M]+ over [M+k]+)
    std::unordered_map<std::string, P1IsoRow, KeyHash> trigger_map;

    // ── Build a set of {gidx||lib_group} pairs where the MONOISOTOPE [N][M]+
    // was found in Phase 1 for that lib_group.  Only these groups are eligible
    // for Phase2_IsoExpansion.  A satellite ([M+k]+_exp, k>=1) from a foreign
    // isotope group — whose monoisotope belongs to a different compound — must
    // NOT trigger expansion; otherwise cross-compound isotope contamination
    // propagates into the output (e.g. [119][M+1]+_exp in Lib_group_4 when
    // [119][M]+ belongs to a feature of Lib_group_7).
    std::unordered_set<std::string, KeyHash> mono_present;  // "{gidx}||{lib_group}"

    // Pass 1: seed from Phase 1 rows that carry an explicit monoisotope iso label.
    // These are features where iso_annotation is non-empty and identifies [M]+ role.
    for (const Rcpp::List& row : rows) {
      std::string iso_ann = Rcpp::as<std::string>(row["iso_annotation"]);
      if (iso_ann.empty()) continue;
      // Monoisotope: has "[M]" with NO "+k" after M inside brackets
      // Format: "[N][M]+" or "[N][M]c+" — detect absence of "[M+"
      size_t p = iso_ann.find(']');
      if (p == std::string::npos) continue;
      size_t q = iso_ann.find('[', p);
      if (q == std::string::npos) continue;
      bool is_mono = !(iso_ann.size() > q + 2 &&
                       iso_ann[q+1] == 'M' &&
                       iso_ann[q+2] == '+');
      if (!is_mono) continue;
      size_t b1 = iso_ann.find('[');
      size_t e1 = iso_ann.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) continue;
      std::string gidx      = iso_ann.substr(b1 + 1, e1 - b1 - 1);
      std::string lib_group = Rcpp::as<std::string>(row["Lib_group"]);
      mono_present.insert(gidx + "||" + lib_group);
    }

    // Pass 2: seed from Phase 1 rows whose iso label was stripped by the multimer
    // stripping guard (iso_annotation empty) but whose mz matches the global
    // iso_group_map monoisotope position for a group that has a _exp satellite in
    // Phase 1 for the same lib_group.
    //
    // Any Phase 1 feature — including Tier2_Unexplained — can serve as a proxy,
    // since all Phase 1 features are already attributed to this compound via
    // library matching. The mz proximity check (within iso_threshold/max_charge
    // + mzabs of the monoisotope) is the primary guard against spurious seeding.
    {
      // First collect all gidx values that need a mono_present seed from _exp
      // satellites already seen, and which don't already have one from pass 1.
      std::unordered_set<std::string, KeyHash> need_seed;  // "{gidx}||{lib_group}"
      for (const Rcpp::List& row : rows) {
        std::string iso_ann = Rcpp::as<std::string>(row["iso_annotation"]);
        if (iso_ann.empty() || iso_ann.find("_exp") == std::string::npos) continue;
        size_t b1 = iso_ann.find('[');
        size_t e1 = iso_ann.find(']');
        if (b1 == std::string::npos || e1 == std::string::npos) continue;
        std::string gidx = iso_ann.substr(b1 + 1, e1 - b1 - 1);
        std::string lib_group = Rcpp::as<std::string>(row["Lib_group"]);
        std::string tkey = gidx + "||" + lib_group;
        if (!mono_present.count(tkey)) need_seed.insert(tkey);
      }

      // For each needed group, find the monoisotope mz from iso_group_map
      // and check whether any Phase 1 row for that lib_group matches it.
      double iso_half = iso_threshold / (double)max_charge;
      double match_tol = iso_half + mzabs;  // conservative: half-spacing + abs floor
      for (const std::string& tkey : need_seed) {
        size_t sep = tkey.find("||");
        if (sep == std::string::npos) continue;
        std::string gidx      = tkey.substr(0, sep);
        std::string lib_group = tkey.substr(sep + 2);

        // Look up monoisotope mz from iso_group_map for this gidx
        auto git = iso_group_map.find(extract_iso_gidx("[" + gidx + "]"));
        if (git == iso_group_map.end()) continue;
        double mono_mz = NA_REAL;
        for (const IsoPoolFeat& ipf : git->second) {
          // Monoisotope label contains "[M]" without "[M+"
          if (ipf.label.find("[M]") != std::string::npos &&
              ipf.label.find("[M+") == std::string::npos) {
            mono_mz = ipf.mz;
            break;
          }
        }
        if (std::isnan(mono_mz)) continue;

        // Check Phase 1 rows: does any row for this lib_group have mz within
        // match_tol of mono_mz?
        for (const Rcpp::List& row : rows) {
          std::string lg = Rcpp::as<std::string>(row["Lib_group"]);
          if (lg != lib_group) continue;
          double rmz = Rcpp::as<double>(row["mz"]);
          if (std::fabs(rmz - mono_mz) <= match_tol) {
            mono_present.insert(tkey);
            break;
          }
        }
      }
    }

    for (const Rcpp::List& row : rows) {
      std::string iso_ann = Rcpp::as<std::string>(row["iso_annotation"]);
      if (iso_ann.empty() || iso_ann.find("_exp") == std::string::npos) continue;

      // Extract gidx from "[N][M...]_exp"
      size_t b1 = iso_ann.find('[');
      size_t e1 = iso_ann.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) continue;
      std::string gidx = iso_ann.substr(b1 + 1, e1 - b1 - 1);

      std::string lib_group = Rcpp::as<std::string>(row["Lib_group"]);

      // ── Compound-scope guard ───────────────────────────────────────────────
      // Only allow expansion if this group's monoisotope was found in Phase 1
      // for the same lib_group.  A satellite from a foreign isotope group
      // (monoisotope absent from this lib_group) is rejected here.
      if (!mono_present.count(gidx + "||" + lib_group)) continue;

      std::string tkey = gidx + "||" + lib_group;

      // Prefer [M]+ as the trigger representative; only replace if not set yet
      bool is_mono = (iso_ann.find("[M]") != std::string::npos &&
                      iso_ann.find("[M+") == std::string::npos);
      auto it = trigger_map.find(tkey);
      if (it == trigger_map.end() || is_mono) {
        double sat_rt = Rcpp::as<double>(row["rt"]);
        // Find best ref_rt for this lib_group: the confirmed_M cluster whose
        // ref_rt is closest to sat_rt (handles multi-cluster same-compound).
        double best_refrt = NA_REAL;
        {
          double best_dist = std::numeric_limits<double>::infinity();
          for (int ci = 0; ci < (int)cm_lib.size(); ++ci) {
            if (cm_lib[ci] != lib_group) continue;
            double crt = cm_ref_rt[ci];
            if (Rcpp::NumericVector::is_na(crt)) continue;
            double dist = std::fabs(crt - (Rcpp::NumericVector::is_na(sat_rt) ? crt : sat_rt));
            if (dist < best_dist) { best_dist = dist; best_refrt = crt; }
          }
        }
        trigger_map[tkey] = {
          Rcpp::as<std::string>(row["feature_name"]),
          Rcpp::as<std::string>(row["Compound_name"]),
          lib_group,
          iso_ann,
          gidx,
          extract_local_idx(iso_ann),  // empty for _exp labels; gidx used as fallback
          Rcpp::as<double>(row["Compound_M"]),
          Rcpp::as<double>(row["Match_peak_score"]),
          Rcpp::as<double>(row["inner_cosine"]),
          Rcpp::as<double>(row["outer_cosine"]),
          Rcpp::as<double>(row["mz"]),
          sat_rt,
          best_refrt
        };
      }
    }

    // ── Seed trigger_map from in-group iso anchors (no _exp) ─────────────
    // Phase 1 features that received an in-group iso label ("[grpN_k][M]+",
    // no "_exp") are monoisotopes confirmed within the library group.
    // Their satellite may exist in the global iso_group_map but was NOT
    // inside the library group — so it never got an _exp label and never
    // entered the trigger loop above.
    //
    // We look up each such anchor in iso_label_map to find its global group
    // index, then seed trigger_map exactly as the _exp path does.
    // The mono_present guard is not needed here: the anchor itself IS in
    // Phase 1 for this lib_group, so compound identity is already confirmed.
    for (const Rcpp::List& row : rows) {
      std::string iso_ann = Rcpp::as<std::string>(row["iso_annotation"]);
      // Only in-group labels: non-empty, no "_exp", contains "[M]" (anchor)
      if (iso_ann.empty()) continue;
      if (iso_ann.find("_exp") != std::string::npos) continue;
      // Must be monoisotope anchor (contains "[M]" but not "[M+")
      bool is_mono = (iso_ann.find("[M]") != std::string::npos &&
                      iso_ann.find("[M+") == std::string::npos);
      if (!is_mono) continue;

      std::string fn_row    = Rcpp::as<std::string>(row["feature_name"]);
      std::string lib_group = Rcpp::as<std::string>(row["Lib_group"]);

      // Look up the global iso group index for this feature in iso_label_map
      auto iit = iso_label_map.find(fn_row);
      if (iit == iso_label_map.end() || iit->second.empty()) continue;
      const std::string& global_lbl = iit->second;
      size_t b1 = global_lbl.find('[');
      size_t e1 = global_lbl.find(']');
      if (b1 == std::string::npos || e1 == std::string::npos) continue;
      std::string gidx = global_lbl.substr(b1 + 1, e1 - b1 - 1);

      // Check whether this global group has any members NOT in phase1_fn_set
      // — only seed if there is actually something to expand.
      auto git = iso_group_map.find(extract_iso_gidx("[" + gidx + "]"));
      if (git == iso_group_map.end()) continue;
      bool has_missing = false;
      for (const IsoPoolFeat& pf : git->second) {
        if (pf.fn != fn_row && !phase1_fn_set.count(pf.fn)) {
          has_missing = true; break;
        }
      }
      if (!has_missing) continue;

      std::string tkey = gidx + "||" + lib_group;
      if (trigger_map.count(tkey)) continue;  // already seeded by _exp path

      double sat_rt = Rcpp::as<double>(row["rt"]);
      double best_refrt = NA_REAL;
      {
        double best_dist = std::numeric_limits<double>::infinity();
        for (int ci = 0; ci < (int)cm_lib.size(); ++ci) {
          if (cm_lib[ci] != lib_group) continue;
          double crt = cm_ref_rt[ci];
          if (Rcpp::NumericVector::is_na(crt)) continue;
          double dist = std::fabs(crt - (Rcpp::NumericVector::is_na(sat_rt) ? crt : sat_rt));
          if (dist < best_dist) { best_dist = dist; best_refrt = crt; }
        }
      }
      // Use the global iso label as the trigger iso_label so Path B label
      // construction can extract the trigger's own k-index correctly.
      trigger_map[tkey] = {
        fn_row,
        Rcpp::as<std::string>(row["Compound_name"]),
        lib_group,
        global_lbl,   // global iso label (no _exp) — gidx already extracted
        gidx,
        extract_local_idx(iso_ann),  // e.g. "66" from "[grp7_66][M]+" — used in emitted labels
        Rcpp::as<double>(row["Compound_M"]),
        Rcpp::as<double>(row["Match_peak_score"]),
        Rcpp::as<double>(row["inner_cosine"]),
        Rcpp::as<double>(row["outer_cosine"]),
        Rcpp::as<double>(row["mz"]),
        sat_rt,
        best_refrt
      };
      // Register monoisotope as present so Path A compound-scope guard passes
      mono_present.insert(tkey);
    }

    // ── Helper lambda: emit one partner row ───────────────────────────────
    // Captures rows, phase1_fn_set, phase1_fn_vec by reference.
    auto emit_partner = [&](const std::string& pfn,
                            double pmz, double prt, double pint,
                            const std::string& plabel,
                            const P1IsoRow& pr) {
      if (phase1_fn_set.insert(pfn).second)
        phase1_fn_vec.push_back(pfn);
      rows.push_back(Rcpp::List::create(
          Rcpp::Named("feature_name")          = pfn,
          Rcpp::Named("mz")                    = pmz,
          Rcpp::Named("rt")                    = prt,
          Rcpp::Named("intensity")             = pint,
          Rcpp::Named("Compound_name")         = pr.compound_name,
          Rcpp::Named("Lib_group")             = pr.lib_group,
          Rcpp::Named("Compound_M")            = pr.M,
          Rcpp::Named("PI_match")              = std::string(""),
          Rcpp::Named("grp_anno_match")      = std::string(""),
          Rcpp::Named("grp_anno_ppm")        = std::string(""),
          Rcpp::Named("grp_anno_mda")        = std::string(""),
          Rcpp::Named("grp_anno_confidence") = std::string(""),
          Rcpp::Named("exp_anno_match")      = std::string(""),
          Rcpp::Named("exp_anno_ppm")        = std::string(""),
          Rcpp::Named("exp_anno_mda")        = std::string(""),
          Rcpp::Named("exp_anno_confidence") = std::string(""),
          Rcpp::Named("NL_match")              = std::string(""),
          Rcpp::Named("iso_annotation")        = plabel,
          Rcpp::Named("reference_tier")       = std::string("Tier1b_exp_Isotope"),
          Rcpp::Named("annotation_phase")      = std::string("Phase2"),
          Rcpp::Named("Match_peak_score")      = pr.match_sc,
          Rcpp::Named("inner_cosine")          = pr.inner_cos,
          Rcpp::Named("outer_cosine")          = pr.outer_cos
      ));
    };

    // ── Expand each {gidx + lib_group} group exactly once ─────────────────
    // emitted tracks {fn||lib_group} pairs to avoid duplicate rows within
    // the same compound context across both paths.
    std::unordered_set<std::string, KeyHash> emitted;

    for (auto& kv : trigger_map) {
      const P1IsoRow& pr = kv.second;

      // ── Path A: iso_group_map lookup ───────────────────────────────────
      // Emits all labeled members of this iso group that are not the
      // triggering feature itself.  We do NOT skip based on phase1_fn_set
      // here: a satellite already annotated under a different compound (e.g.
      // it matched a different library group) must still appear under this
      // compound too — the same feature appearing under multiple compounds is
      // valid and consistent with how adduct matches work.
      // Deduplication within the same {fn + lib_group} is handled by the
      // emitted set so we never emit the same row twice.
      std::unordered_set<std::string, KeyHash> path_a_fns;
      {
        auto git = iso_group_map.find(extract_iso_gidx("[" + pr.gidx + "]"));
        if (git != iso_group_map.end()) {

          // ── Compound-scoped isotope group check ───────────────────────────
          // Check_Iso_Charge assigns GLOBAL group indices across the full pool.
          // A trigger feature with label "[119][M+1]+_exp" belongs to group 119
          // which may have been formed from features of OTHER compounds.
          // Emitting all group 119 members would pull foreign-compound isotope
          // satellites into this compound's annotation.
          //
          // Guard: the monoisotope [N][M]+ of this group must itself be in
          // phase1_fn_set (already annotated for this compound) OR must be the
          // trigger feature itself (pr.fn).  If the monoisotope is absent from
          // this compound's annotation context, the group is foreign — skip it.
          std::string mono_fn;
          for (const IsoPoolFeat& pf : git->second) {
            // Monoisotope has label "[N][M]+" — no "+k" inside brackets
            size_t p = pf.label.find(']');
            if (p == std::string::npos) continue;
            size_t q = pf.label.find('[', p);
            if (q == std::string::npos) continue;
            bool is_mono = !(pf.label.size() > q + 2 &&
                             pf.label[q+1] == 'M' &&
                             pf.label[q+2] == '+');
            if (is_mono) { mono_fn = pf.fn; break; }
          }

          // If monoisotope is not the trigger itself and not in phase1_fn_set,
          // this iso group is foreign to this compound — skip entirely.
          bool mono_ok = mono_fn.empty()          // no monoisotope found (odd edge case — allow)
                      || mono_fn == pr.fn
                      || phase1_fn_set.count(mono_fn);
          if (!mono_ok) goto path_b;  // skip Path A for this trigger

          for (const IsoPoolFeat& pf : git->second) {
            if (pf.fn == pr.fn) continue;   // skip self (the triggering feature)

            // Skip only if this feature already has a Phase 1 iso label for
            // this same lib_group — it has already been correctly handled.
            // Features that are Tier2/Tier1d/Tier1c in Phase 1 with no iso
            // label must NOT be blocked: they are unexplained and legitimately
            // need an isotope satellite label from Path A.
            // Cross-compound features (same fn, different lib_group) are also
            // allowed through — multi-compound assignment is by design.
            {
              bool already_has_iso = false;
              for (const Rcpp::List& r : rows) {
                if (Rcpp::as<std::string>(r["feature_name"]) == pf.fn &&
                    Rcpp::as<std::string>(r["Lib_group"])     == pr.lib_group) {
                  std::string iso = Rcpp::as<std::string>(r["iso_annotation"]);
                  if (!iso.empty()) { already_has_iso = true; }
                  break;
                }
              }
              if (already_has_iso) continue;
            }

            std::string key = pf.fn + "||" + pr.lib_group;
            if (!emitted.insert(key).second) continue;
            path_a_fns.insert(pf.fn);

            // Add to phase1_fn_set only if not already there, so Phase 2
            // M-expansion does not re-emit it as a fresh adduct candidate.
            if (phase1_fn_set.insert(pf.fn).second)
              phase1_fn_vec.push_back(pf.fn);

            // ── Path A label construction ──────────────────────────────────
            // Normally: prefix pf.label with lib group index.
            // Special case: if pf.label marks a MONOISOTOPE ([M]c+ with c>=2)
            // but pf.mz is within satellite spacing of the trigger anchor (pr.mz),
            // the global run assigned it a new group because the anchor was labeled
            // singly charged — rewrite to use the anchor's group index and the
            // correct [M+k]c+ satellite position instead.
            std::string pa_lbl = pf.label;
            {
              // Check if pf is a global-run monoisotope at multiply-charged
              // satellite spacing from the trigger anchor.
              bool rewritten = false;
              if (!Rcpp::NumericVector::is_na(pr.mz) &&
                  !Rcpp::NumericVector::is_na(pf.mz)) {
                // Extract ion type from pf.label: second bracket pair content
                size_t e1p = pf.label.find(']');
                size_t b2p = (e1p != std::string::npos) ?
                             pf.label.find('[', e1p + 1) : std::string::npos;
                if (b2p != std::string::npos) {
                  std::string ion_part = pf.label.substr(b2p + 1);
                  // Is it a monoisotope at c>=2? e.g. "M]2+" or "M]3+"
                  bool is_mc_mono = (ion_part.size() >= 2 &&
                                     ion_part[0] == 'M' &&
                                     ion_part[1] != '+');
                  if (is_mc_mono) {
                    double mz_diff = pf.mz - pr.mz;
                    for (int cc = 2; cc <= max_charge && !rewritten; ++cc) {
                      for (int mm = 1; mm <= max_isotope && !rewritten; ++mm) {
                        double expected = mm * iso_threshold / (double)cc;
                        if (std::fabs(mz_diff - expected) <= mzabs * 3.0) {
                          // Build corrected label: [gidx][M+mm]cc+_exp
                          std::string grp_num = pr.lib_group;
                          size_t us = grp_num.rfind('_');
                          if (us != std::string::npos)
                            grp_num = grp_num.substr(us + 1);
                          std::string c_str = (cc > 1) ? std::to_string(cc) : "";
                          // pr.gidx may be "grp3_32" (in-group format) —
                          // strip the grpN_ prefix to get just the numeric index.
                          // Prefer pr.local_idx (in-group index, e.g. "32") when set.
                          std::string clean_gidx = pr.local_idx.empty()
                            ? extract_iso_gidx("[" + pr.gidx + "]")
                            : pr.local_idx;
                          if (clean_gidx.empty()) clean_gidx = pr.gidx;
                          pa_lbl = "[grp" + grp_num + "_" + clean_gidx + "]" +
                                   "[M+" + std::to_string(mm) + "]" +
                                   c_str + "+";
                          rewritten = true;
                        }
                      }
                    }
                  }
                }
              }
              if (!rewritten) {
                // Standard label: prefix with lib group, using the trigger's
                // local in-group index when available (local_idx), otherwise
                // the global gidx from pf.label.  This ensures emitted satellite
                // labels share the same index as the anchor's Phase 1 iso_annotation.
                size_t b1 = pa_lbl.find('['), e1 = pa_lbl.find(']');
                if (b1 != std::string::npos && e1 != std::string::npos) {
                  std::string grp_num = pr.lib_group;
                  size_t us = grp_num.rfind('_');
                  if (us != std::string::npos) grp_num = grp_num.substr(us + 1);
                  std::string use_idx = pr.local_idx.empty()
                    ? pa_lbl.substr(b1 + 1, e1 - b1 - 1)  // global gidx from pf.label
                    : pr.local_idx;                          // local in-group index
                  pa_lbl = "[grp" + grp_num + "_" + use_idx + "]" +
                           pa_lbl.substr(e1 + 1);
                }
              }
            }
            emit_partner(pf.fn, pf.mz, pf.rt, pf.intensity,
                         pa_lbl + "_exp", pr);
          }
        }
      }
      path_b:;

      // ── Path B: direct pool scan for unlabeled partners ────────────────
      // Covers satellites that Check_Iso_Charge could not label because the
      // m=1 intensity check failed (adduct form intensity ≠ true monoisotope).
      // Scans pool in BOTH directions from the triggering feature's mz:
      //   forward  (mz > trigger): looks for [M+k]+ satellites
      //   backward (mz < trigger): looks for lower-m members or monoisotope
      // Uses a 30 mDa window (3 × mzabs) — wider than the iso run — because
      // we already have strong compound identity evidence from the adduct match.
      if (!Rcpp::NumericVector::is_na(pr.mz) && n_pool > 0) {
        for (int p = 0; p < n_pool; ++p) {
          const std::string& fn_p = pool_fn_vec[p];
          if (fn_p == pr.fn) continue;
          // For phase1 features: only skip if they already have an iso label —
          // UNLESS Path B determines this is a multiply-charged satellite whose
          // global-fallback label used the wrong iso group (a different gidx than
          // the current trigger's group).  In that case we fall through and allow
          // the group-override branch (c >= 2 && pool_is_mono) to fire, then
          // update the Phase 1 row's iso_annotation in place rather than emitting
          // a new row.  has_iso_row_idx records the rows[] index for that update.
          int has_iso_row_idx = -1;
          if (phase1_fn_set.count(fn_p)) {
            bool has_iso = false;
            for (int ri = 0; ri < (int)rows.size(); ++ri) {
              if (Rcpp::as<std::string>(rows[ri]["feature_name"]) == fn_p &&
                  Rcpp::as<std::string>(rows[ri]["Lib_group"])     == pr.lib_group) {
                std::string iso = Rcpp::as<std::string>(rows[ri]["iso_annotation"]);
                if (!iso.empty()) { has_iso = true; has_iso_row_idx = ri; break; }
              }
            }
            if (has_iso && has_iso_row_idx < 0) continue; // has iso under different lib_group — skip
            if (has_iso && has_iso_row_idx >= 0) {
              // Check whether the existing label belongs to the same iso group as
              // the current trigger.  If the groups already match, no correction
              // needed — skip.  If they differ, fall through so the group-override
              // branch can correct it.
              const std::string& existing_iso =
                Rcpp::as<std::string>(rows[has_iso_row_idx]["iso_annotation"]);
              std::string norm_gidx_trigger = extract_iso_gidx("[" + pr.gidx + "]");
              if (norm_gidx_trigger.empty()) norm_gidx_trigger = pr.gidx;
              size_t eb1 = existing_iso.find('['), ee1 = existing_iso.find(']');
              std::string existing_gidx = (eb1 != std::string::npos && ee1 != std::string::npos)
                ? extract_iso_gidx(existing_iso.substr(eb1, ee1 - eb1 + 1))
                : "";
              if (existing_gidx == norm_gidx_trigger) continue; // same group — no correction needed
              // Different group: fall through to allow group-override correction
            }
          }
          if (path_a_fns.count(fn_p)) continue;
          if (Rcpp::NumericVector::is_na(pool_mz_vec[p])) continue;

          // RT gate: use the [M+H]+ anchor RT (pr.ref_rt) to correctly scope
          // the pool scan to features co-eluting with the adduct-matched
          // precursor — not the satellite's own (potentially shifted) RT.
          // Fall back to pr.rt when ref_rt is NA (no confirmed_M entry).
          double rt_anchor = Rcpp::NumericVector::is_na(pr.ref_rt) ? pr.rt : pr.ref_rt;
          if (!Rcpp::NumericVector::is_na(rt_anchor) &&
              !Rcpp::NumericVector::is_na(pool_rt_vec[p])) {
              if (std::fabs(pool_rt_vec[p] - rt_anchor) > rt_exp_threshold) continue;
          }

          double mz_diff = pool_mz_vec[p] - pr.mz;  // can be negative (lower member)

          // Test all valid {m, c} combinations in both directions
          bool matched = false;
          std::string sat_label;
          for (int c = 1; c <= max_charge && !matched; ++c) {
            for (int m = -max_isotope; m <= max_isotope && !matched; ++m) {
              if (m == 0) continue;
              double expected = m * iso_threshold / (double)c;
              if (std::fabs(mz_diff - expected) <= mzabs * 3.0) {
                matched = true;
                // Build label: m > 0 → higher satellite, m < 0 → lower member
                int abs_m = std::abs(m);
                std::ostringstream lbl;
                // Use the pool feature's own global iso group index if available.
                // Inheriting pr.gidx (the trigger's group index) is wrong when
                // the trigger belongs to a different iso group than the pool
                // feature — e.g. trigger is group 343 (FT2279) but pool feature
                // FT1771 belongs to group 258 (paired with FT1768).
                // If the pool feature belongs to a DIFFERENT group, skip it
                // entirely — coincidental mz spacing, not a real satellite.
                // Fall back to pr.gidx only when pool feature has no global label.
                {
                  // Normalise pr.gidx: strip grpN_ prefix if present
                  std::string norm_gidx = extract_iso_gidx("[" + pr.gidx + "]");
                  if (norm_gidx.empty()) norm_gidx = pr.gidx;
                  std::string use_gidx = norm_gidx;
                  auto pool_iit = iso_label_map.find(fn_p);
                  if (pool_iit != iso_label_map.end() && !pool_iit->second.empty()) {
                    size_t pb1 = pool_iit->second.find('[');
                    size_t pe1 = pool_iit->second.find(']');
                    if (pb1 != std::string::npos && pe1 != std::string::npos) {
                      std::string pool_gidx = pool_iit->second.substr(pb1 + 1, pe1 - pb1 - 1);
                      if (pool_gidx != norm_gidx) {
                        // Groups differ — could be a multiply-charged satellite that
                        // Check_Iso_Charge missed because the anchor was labeled as
                        // singly charged (c=1) while the satellite is c>=2.
                        // Allow if:
                        //   (a) c >= 2 (fractional Da spacing, e.g. +0.5 for c=2)
                        //   (b) pool feature is labeled as a MONOISOTOPE in its own
                        //       global group (not a satellite of something else),
                        //       meaning the global run assigned it a new group rather
                        //       than recognising it as a satellite of the anchor.
                        // In this case adopt pr.gidx and build the correct label.
                        bool pool_is_mono = false;
                        {
                          const std::string& pl = pool_iit->second;
                          size_t e1p = pl.find(']');
                          size_t b2p = (e1p != std::string::npos) ?
                                       pl.find('[', e1p+1) : std::string::npos;
                          if (b2p != std::string::npos) {
                            std::string ion_part = pl.substr(b2p + 1);
                            // Mono: starts "M]" or "M]c" (no "M+" in ion type)
                            pool_is_mono = (ion_part.size() >= 2 &&
                                           ion_part[0] == 'M' &&
                                           ion_part[1] != '+');
                          }
                        }
                        if (c >= 2 && pool_is_mono) {
                          // Override: adopt anchor group index (normalised), treat
                          // pool feature as a multiply-charged satellite of anchor.
                          use_gidx = norm_gidx;
                        } else {
                          matched = false; break;
                        }
                      } else {
                        use_gidx = pool_gidx;
                      }
                    }
                  }
                  lbl << "[" << (pr.local_idx.empty() ? use_gidx : pr.local_idx) << "]";
                }
                if (m > 0)
                  lbl << "[M+" << abs_m << "]";
                else
                  // This pool feature is m steps BELOW the trigger:
                  // e.g. trigger is [M+2]+, pool feature is [M]+ (m = -2)
                  // We compute its actual isotope index relative to the group.
                  // The trigger's own isotope index within the group is
                  // extracted from its label "[N][M+k]+_exp" → k.
                  // The pool feature's index = trigger_k + m.
                {
                  // Extract trigger's own m-index from iso_label e.g. "[N][M+2]+_exp"
                  int trigger_k = 0;
                  {
                    std::string& lbl_ref = const_cast<std::string&>(pr.iso_label);
                    size_t plus_pos = lbl_ref.find("[M+");
                    if (plus_pos != std::string::npos) {
                      size_t num_start = plus_pos + 3;
                      size_t num_end   = lbl_ref.find(']', num_start);
                      if (num_end != std::string::npos)
                        trigger_k = std::stoi(lbl_ref.substr(num_start,
                                                             num_end - num_start));
                    }
                    // "[N][M]+" → trigger_k stays 0 (monoisotope)
                  }
                  int pool_k = trigger_k + m;  // m is negative here
                  if (pool_k < 0) { matched = false; break; }  // below monoisotope — invalid
                  if (pool_k == 0)
                    lbl << "[M]";
                  else
                    lbl << "[M+" << pool_k << "]";
                }
                if (c > 1) lbl << c;
                lbl << "+_exp";
                // Prefix with lib group number for traceability
                {
                  std::string raw = lbl.str();
                  size_t b1 = raw.find('['), e1 = raw.find(']');
                  if (b1 != std::string::npos && e1 != std::string::npos) {
                    std::string grp_num = pr.lib_group;
                    size_t us = grp_num.rfind('_');
                    if (us != std::string::npos) grp_num = grp_num.substr(us + 1);
                    raw = "[grp" + grp_num + "_" +
                          raw.substr(b1 + 1, e1 - b1 - 1) + "]" +
                          raw.substr(e1 + 1);
                  }
                  sat_label = raw;
                }
              }
            }
          }
          if (!matched) continue;

          // ── Sequential isotope check ──────────────────────────────────────
          // Extract the target isotope index k from sat_label "[N][M+k]+_exp".
          // For k > 1, all intermediate steps [M+1]+..[M+k-1]+ must already
          // be present in phase1_fn_set or emitted — a gap (e.g. [M+4]+ with
          // no [M+2]+, [M+3]+) is physically implausible and rejected.
          {
            int target_k = 0;
            {
              size_t p1 = sat_label.find("[M+");
              if (p1 != std::string::npos) {
                size_t n0 = p1 + 3;
                size_t n1 = sat_label.find(']', n0);
                if (n1 != std::string::npos)
                  target_k = std::stoi(sat_label.substr(n0, n1 - n0));
              }
              // "[N][M]+_exp" → target_k stays 0 (monoisotope, always allowed)
            }
            if (target_k > 1) {
              // Check that every intermediate [M+1]+..[M+k-1]+ was already
              // emitted for this iso group.  We scan emitted set for each step.
              // The emitted key format is "fn||lib_group"; we need to check
              // whether ANY fn in the group has label "[N][M+j]+" for each j.
              // Simpler proxy: check iso_group_map for the group — if the group
              // has fewer than (target_k) labeled members (excluding monoisotope),
              // the intermediate steps are missing.
              bool gap_found = false;
              auto git = iso_group_map.find(extract_iso_gidx("[" + pr.gidx + "]"));
              if (git == iso_group_map.end()) {
                gap_found = true;  // group unknown — can't verify continuity
              } else {
                // Build set of already-known isotope indices for this group
                std::unordered_set<int> known_k;
                for (const IsoPoolFeat& pf : git->second) {
                  // Extract k from pf.label
                  size_t pp = pf.label.find("[M+");
                  if (pp == std::string::npos) { known_k.insert(0); continue; }
                  size_t n0 = pp + 3;
                  size_t n1 = pf.label.find(']', n0);
                  if (n1 != std::string::npos)
                    known_k.insert(std::stoi(pf.label.substr(n0, n1 - n0)));
                }
                // Also check already-emitted Path B features for this group
                // by scanning phase1_fn_set labels — approximation: accept
                // if known_k contains all 1..(target_k-1)
                for (int j = 1; j < target_k; ++j) {
                  if (!known_k.count(j)) { gap_found = true; break; }
                }
              }
              if (gap_found) continue;  // reject: intermediate peak missing
            }
          }

          std::string key = fn_p + "||" + pr.lib_group;
          if (!emitted.insert(key).second) continue;

          // If this feature was already emitted in Phase 1 with a wrong-group
          // iso label (has_iso_row_idx >= 0), update that row in place instead
          // of emitting a new one.  This corrects labels like [grp3_538][M]2+_exp
          // → [grp3_537][M+1]2+_exp when the global run assigned the feature to
          // its own group (538) but Path B identifies it as a satellite of the
          // current trigger's group (537).
          if (has_iso_row_idx >= 0) {
            rows[has_iso_row_idx]["iso_annotation"] = sat_label;
            // Upgrade tier if still Tier2_Unexplained
            std::string cur_tier = Rcpp::as<std::string>(rows[has_iso_row_idx]["reference_tier"]);
            if (cur_tier == "Tier2_Unexplained")
              rows[has_iso_row_idx]["reference_tier"] = std::string("Tier1b_exp_Isotope");
          } else {
            emit_partner(fn_p, pool_mz_vec[p], pool_rt_vec[p], pool_int_vec[p],
                         sat_label, pr);
          }
        }
      }
    }

    // ── Path C: direct pool scan for Phase 1 features with NO iso label ───────
    // These features were missed by the global Check_Iso_Charge run (consumed
    // as someone else's satellite internally, intensity ratio failed, or simply
    // not paired) — so no _exp label propagated into Phase 1 rows, trigger_map
    // has no entry for them, and Paths A/B never fire.
    //
    // IMPORTANT: Only Tier1a_AnnoMatch features are valid anchors here.
    // Tier1b/1c/1d features are isotope satellites, MS2 fragments, or NL
    // features — not confirmed monoisotopes. Using them as Path C anchors
    // produces spurious satellite labels from coincidental 1 Da spacings.
    // Tier1a features have a confirmed adduct match and are genuine candidates
    // for being the monoisotope of a real isotope group.
    if (n_pool > 0) {
      std::unordered_set<std::string, KeyHash> trigger_fns;
      for (auto& kv : trigger_map) trigger_fns.insert(kv.second.fn);

      for (const Rcpp::List& row : rows) {
        std::string iso_ann = Rcpp::as<std::string>(row["iso_annotation"]);
        // Admit Tier1a anchors that either have no iso label, OR have an
        // in-group monoisotope label ("[M]+" role, no "[M+").  Satellites and
        // _exp-labeled features are excluded — only confirmed adduct-matched
        // monoisotopes are valid Path C anchors.
        if (!iso_ann.empty()) {
          // Allow only if it is a monoisotope (no "[M+" after the group bracket)
          // and NOT an _exp label (those are handled by Paths A/B).
          if (iso_ann.find("_exp") != std::string::npos) continue;
          // Must be monoisotope: second bracket pair starts "M]", not "M+"
          size_t e1c = iso_ann.find(']');
          size_t b2c = (e1c != std::string::npos) ? iso_ann.find('[', e1c+1) : std::string::npos;
          if (b2c == std::string::npos) continue;
          std::string ion_part_c = iso_ann.substr(b2c + 1);
          bool is_mono_anchor = (ion_part_c.size() >= 2 &&
                                 ion_part_c[0] == 'M' && ion_part_c[1] != '+');
          if (!is_mono_anchor) continue;
        }
        std::string ann_phase = Rcpp::as<std::string>(row["annotation_phase"]);
        if (ann_phase != "Phase1") continue;
        // Only Tier1a_AnnoMatch features are valid isotope monoisotope anchors.
        std::string tier = Rcpp::as<std::string>(row["reference_tier"]);
        if (tier.find("Tier1a_AnnoMatch") == std::string::npos) continue;
        std::string fn_anchor = Rcpp::as<std::string>(row["feature_name"]);
        if (!phase1_fn_set.count(fn_anchor)) continue;
        if (trigger_fns.count(fn_anchor)) continue;

        std::string lib_group     = Rcpp::as<std::string>(row["Lib_group"]);
        std::string compound_name = Rcpp::as<std::string>(row["Compound_name"]);
        double anchor_mz = Rcpp::as<double>(row["mz"]);
        double anchor_rt = Rcpp::as<double>(row["rt"]);
        double anchor_M  = Rcpp::as<double>(row["Compound_M"]);
        double anchor_sc = Rcpp::as<double>(row["Match_peak_score"]);
        double anchor_ic = Rcpp::as<double>(row["inner_cosine"]);
        double anchor_oc = Rcpp::as<double>(row["outer_cosine"]);
        if (Rcpp::NumericVector::is_na(anchor_mz)) continue;

        double ref_rt = NA_REAL;
        {
          double best_dist = std::numeric_limits<double>::infinity();
          for (int ci = 0; ci < (int)cm_lib.size(); ++ci) {
            if (cm_lib[ci] != lib_group) continue;
            double crt = cm_ref_rt[ci];
            if (Rcpp::NumericVector::is_na(crt)) continue;
            double dist = std::fabs(crt - (Rcpp::NumericVector::is_na(anchor_rt) ? crt : anchor_rt));
            if (dist < best_dist) { best_dist = dist; ref_rt = crt; }
          }
        }
        double rt_gate = Rcpp::NumericVector::is_na(ref_rt) ? anchor_rt : ref_rt;

        std::string grp_num = lib_group;
        size_t us = grp_num.rfind('_');
        if (us != std::string::npos) grp_num = grp_num.substr(us + 1);

        // Extract the iso group index for label construction from the anchor's
        // own iso_annotation when available (e.g. "[grp3_32][M]+" → "32").
        // This is more accurate than the old mz_idx fallback (rounded mz) and
        // produces labels consistent with the in-group Phase 1 iso labels.
        // Fall back to rounded mz only if the anchor has no iso label.
        std::string anchor_iso_idx;
        {
          std::string anchor_iso_ann = Rcpp::as<std::string>(row["iso_annotation"]);
          if (!anchor_iso_ann.empty()) {
            // Format: "[grpN_IDX][M]..." — extract IDX
            size_t ai_b1 = anchor_iso_ann.find('[');
            size_t ai_e1 = anchor_iso_ann.find(']');
            if (ai_b1 != std::string::npos && ai_e1 != std::string::npos) {
              std::string inner = anchor_iso_ann.substr(ai_b1 + 1, ai_e1 - ai_b1 - 1);
              // Strip grpN_ prefix
              size_t ai_us = inner.rfind('_');
              if (ai_us != std::string::npos) {
                std::string prefix = inner.substr(0, ai_us);
                if (prefix.size() >= 3 &&
                    prefix[0]=='g' && prefix[1]=='r' && prefix[2]=='p')
                  anchor_iso_idx = inner.substr(ai_us + 1);
              }
              if (anchor_iso_idx.empty()) anchor_iso_idx = inner;
            }
          }
          if (anchor_iso_idx.empty()) {
            // Fallback: use anchor's global iso label index
            auto anch_iit = iso_label_map.find(fn_anchor);
            if (anch_iit != iso_label_map.end() && !anch_iit->second.empty()) {
              anchor_iso_idx = extract_iso_gidx(anch_iit->second);
            }
          }
          if (anchor_iso_idx.empty()) {
            // Last resort: rounded mz (original behaviour)
            anchor_iso_idx = std::to_string(static_cast<int>(std::round(anchor_mz)));
          }
        }

        // anchor's global iso group index (for pool feature group-match guard)
        std::string anchor_global_gidx;
        {
          auto anch_iit = iso_label_map.find(fn_anchor);
          if (anch_iit != iso_label_map.end() && !anch_iit->second.empty())
            anchor_global_gidx = extract_iso_gidx(anch_iit->second);
        }

        for (int p = 0; p < n_pool; ++p) {
          const std::string& fn_p = pool_fn_vec[p];
          if (fn_p == fn_anchor) continue;
          if (Rcpp::NumericVector::is_na(pool_mz_vec[p])) continue;

          // For Phase 1 features with an iso label: allow in-place correction
          // when their existing label belongs to a different iso group than the
          // anchor (same logic as Path B).  has_iso_row_idx_c records the row.
          int has_iso_row_idx_c = -1;
          if (phase1_fn_set.count(fn_p)) {
            bool has_iso_c = false;
            for (int ri = 0; ri < (int)rows.size(); ++ri) {
              if (Rcpp::as<std::string>(rows[ri]["feature_name"]) == fn_p &&
                  Rcpp::as<std::string>(rows[ri]["Lib_group"])     == lib_group) {
                std::string iso_c = Rcpp::as<std::string>(rows[ri]["iso_annotation"]);
                if (!iso_c.empty()) { has_iso_c = true; has_iso_row_idx_c = ri; break; }
              }
            }
            if (has_iso_c && has_iso_row_idx_c < 0) continue; // iso under different lib_group
            if (has_iso_c && has_iso_row_idx_c >= 0) {
              // Same-group check: if existing label already matches anchor's group, skip
              const std::string& existing_iso_c =
                Rcpp::as<std::string>(rows[has_iso_row_idx_c]["iso_annotation"]);
              size_t ec_b1 = existing_iso_c.find('['), ec_e1 = existing_iso_c.find(']');
              std::string existing_gidx_c = (ec_b1 != std::string::npos && ec_e1 != std::string::npos)
                ? extract_iso_gidx(existing_iso_c.substr(ec_b1, ec_e1 - ec_b1 + 1))
                : "";
              if (existing_gidx_c == anchor_iso_idx) continue; // already correct
              // Different group: fall through for in-place correction
            } else if (!has_iso_c) {
              continue; // in phase1_fn_set but no iso label for this lib_group — skip
            }
          }

          if (!Rcpp::NumericVector::is_na(rt_gate) &&
              !Rcpp::NumericVector::is_na(pool_rt_vec[p])) {
            if (std::fabs(pool_rt_vec[p] - rt_gate) > rt_exp_threshold) continue;
          }

          double mz_diff = pool_mz_vec[p] - anchor_mz;
          if (mz_diff <= 0) continue;

          bool matched = false;
          std::string sat_label;
          // For in-place correction of a Phase 1 feature, we need the group-
          // override check: pool feature must be a monoisotope in its own global
          // group (mis-assigned by the global run) at c >= 2 spacing.
          bool is_inplace_correction = (has_iso_row_idx_c >= 0);
          for (int c = 1; c <= max_charge && !matched; ++c) {
            for (int m = 1; m <= max_isotope && !matched; ++m) {
              double expected = m * iso_threshold / (double)c;
              if (std::fabs(mz_diff - expected) <= mzabs * 1.5) {
                // For in-place correction, require c >= 2 and pool feature must
                // be a monoisotope in its own global group (same logic as Path B).
                if (is_inplace_correction) {
                  if (c < 2) continue;
                  auto pool_iit = iso_label_map.find(fn_p);
                  if (pool_iit == iso_label_map.end() || pool_iit->second.empty()) continue;
                  const std::string& pl = pool_iit->second;
                  size_t e1p = pl.find(']');
                  size_t b2p = (e1p != std::string::npos) ? pl.find('[', e1p+1) : std::string::npos;
                  if (b2p == std::string::npos) continue;
                  std::string ion_part = pl.substr(b2p + 1);
                  bool pool_is_mono_c = (ion_part.size() >= 2 &&
                                         ion_part[0] == 'M' && ion_part[1] != '+');
                  if (!pool_is_mono_c) continue;
                }
                matched = true;
                std::ostringstream lbl;
                lbl << "[grp" << grp_num << "_" << anchor_iso_idx << "][M+" << m << "]";
                if (c > 1) lbl << c;
                lbl << "+_exp";
                sat_label = lbl.str();
              }
            }
          }
          if (!matched) continue;

          // Guard: skip if pool feature's global iso label belongs to a DIFFERENT
          // group than the anchor's, UNLESS in-place correction already approved above.
          if (!is_inplace_correction && !anchor_global_gidx.empty()) {
            auto pool_iit = iso_label_map.find(fn_p);
            if (pool_iit != iso_label_map.end() && !pool_iit->second.empty()) {
              std::string pool_gidx = extract_iso_gidx(pool_iit->second);
              if (pool_gidx != anchor_global_gidx) continue;
            }
          }

          std::string key = fn_p + "||" + lib_group;
          if (!emitted.insert(key).second) continue;

          if (is_inplace_correction) {
            // Update the existing Phase 1 row in place
            rows[has_iso_row_idx_c]["iso_annotation"] = sat_label;
            std::string cur_tier_c = Rcpp::as<std::string>(rows[has_iso_row_idx_c]["reference_tier"]);
            if (cur_tier_c == "Tier2_Unexplained")
              rows[has_iso_row_idx_c]["reference_tier"] = std::string("Tier1b_exp_Isotope");
          } else {
            if (phase1_fn_set.insert(fn_p).second)
              phase1_fn_vec.push_back(fn_p);
            rows.push_back(Rcpp::List::create(
              Rcpp::Named("feature_name")        = fn_p,
              Rcpp::Named("mz")                  = pool_mz_vec[p],
              Rcpp::Named("rt")                  = pool_rt_vec[p],
              Rcpp::Named("intensity")           = pool_int_vec[p],
              Rcpp::Named("Compound_name")       = compound_name,
              Rcpp::Named("Lib_group")           = lib_group,
              Rcpp::Named("Compound_M")          = anchor_M,
              Rcpp::Named("PI_match")            = std::string(""),
              Rcpp::Named("grp_anno_match")      = std::string(""),
              Rcpp::Named("grp_anno_ppm")        = std::string(""),
              Rcpp::Named("grp_anno_mda")        = std::string(""),
              Rcpp::Named("grp_anno_confidence") = std::string(""),
              Rcpp::Named("exp_anno_match")      = std::string(""),
              Rcpp::Named("exp_anno_ppm")        = std::string(""),
              Rcpp::Named("exp_anno_mda")        = std::string(""),
              Rcpp::Named("exp_anno_confidence") = std::string(""),
              Rcpp::Named("NL_match")            = std::string(""),
              Rcpp::Named("iso_annotation")      = sat_label,
              Rcpp::Named("reference_tier")      = std::string("Tier1b_exp_Isotope"),
              Rcpp::Named("annotation_phase")    = std::string("Phase2"),
              Rcpp::Named("Match_peak_score")    = anchor_sc,
              Rcpp::Named("inner_cosine")        = anchor_ic,
              Rcpp::Named("outer_cosine")        = anchor_oc
            ));
          }
        }
      }
    }
  }

  // ── build iso_labels data.frame for return ────────────────────────────────
  std::vector<std::string> il_fn;
  std::vector<std::string> il_label;
  il_fn.reserve(iso_label_map.size());
  il_label.reserve(iso_label_map.size());
  for (auto& kv : iso_label_map) {
    il_fn.push_back(kv.first);
    il_label.push_back(kv.second);
  }

  // ── assemble return value ─────────────────────────────────────────────────
  Rcpp::List rows_out(rows.size());
  for (int i = 0; i < (int)rows.size(); ++i) rows_out[i] = rows[i];

  Rcpp::DataFrame confirmed_M_df = Rcpp::DataFrame::create(
    Rcpp::Named("M")             = cm_M,
    Rcpp::Named("compound_name") = cm_name,
    Rcpp::Named("lib_group")     = cm_lib,
    Rcpp::Named("ref_rt")        = cm_ref_rt,
    Rcpp::Named("stringsAsFactors") = false
  );

  Rcpp::DataFrame iso_labels_df = Rcpp::DataFrame::create(
    Rcpp::Named("feature_name") = il_fn,
    Rcpp::Named("iso_label")    = il_label,
    Rcpp::Named("stringsAsFactors") = false
  );

  Rcpp::CharacterVector phase1_fn_r(phase1_fn_vec.begin(), phase1_fn_vec.end());

  return Rcpp::List::create(
    Rcpp::Named("rows")             = rows_out,
    Rcpp::Named("confirmed_M")      = confirmed_M_df,
    Rcpp::Named("iso_labels")       = iso_labels_df,
    Rcpp::Named("phase1_features")  = phase1_fn_r
  );
}


// ═══════════════════════════════════════════════════════════════════════════
// EXPORTED FUNCTION 2 — Phase 2: M-expansion on the full feature pool
//
// All Phase 2 features are discovered via adduct match in the RT-window pool,
// so all evidence here is "ext" (extended/expansion) scope by definition.
//
// Confidence tier ladder for Phase 2 (priority order):
//   Tier1a_exp_AnnoMatch    adduct match in RT-window pool  [priority 6]
//   Tier1b_exp_Isotope        isotope partner in global pool  [priority 7]
//   Tier1c_MS2Fragment_exp    MS2 fragment match (no Phase 1 group) [priority 8]
//   Tier1d_exp_NeutralLoss    NL within expansion group       [priority 9]
//
// IMPORTANT: Tier1d_exp_NeutralLoss fires ONLY when no adduct match exists.
// NL is supplementary evidence stored in NL_match; it never overwrites an
// adduct-confirmed tier (Tier1a_exp_AnnoMatch is always stronger).
//
// Inputs:
//   FT        — full feature pool table
//   confirmed_M_df — data.frame from Phase 1 ($confirmed_M):
//                    columns M, compound_name, lib_group, ref_rt
//   iso_labels_df  — data.frame from Phase 1 ($iso_labels):
//                    columns feature_name, iso_label
//                    (global iso labels already computed; not re-run here)
//   phase1_features — CharacterVector of feature names already annotated in
//                     Phase 1; these are skipped in the pool scan
//   adducts_df, neutral_loss_table, polarity — same as Phase 1
//   mz_threshold   — absolute Da floor (same as Phase 1)
//   ppm_exp        — ppm gate for Phase 2 adduct search (default 5)
//   nl_ppm         — ppm gate for NL within expansion group (default 5)
//   rt_exp_threshold — RT window around ref_rt (default 5 s)
//
// Returns a List of Phase 2 annotation rows (same column schema as Phase 1).
// ═══════════════════════════════════════════════════════════════════════════

// [[Rcpp::export]]
Rcpp::List annotate_librarymatch_phase2_cpp(
    Rcpp::DataFrame FT,
    Rcpp::DataFrame confirmed_M_df,
    Rcpp::DataFrame iso_labels_df,
    Rcpp::CharacterVector phase1_features,
    Rcpp::DataFrame adducts_df,
    Rcpp::DataFrame neutral_loss_table,
    std::string     polarity,
    double mz_threshold     = 0.005,
    double ppm_exp          = 5.0,
    double nl_ppm           = 5.0,
    double rt_exp_threshold = 5.0)
{
  // ── normalise polarity ────────────────────────────────────────────────────
  polarity = tolower_copy(polarity);
  if      (polarity == "positive" || polarity == "+" || polarity == "pos") polarity = "pos";
  else if (polarity == "negative" || polarity == "-" || polarity == "neg") polarity = "neg";

  // ── adduct table ──────────────────────────────────────────────────────────
  std::vector<std::string> filtered_ions, filtered_formulas;
  build_adduct_table(adducts_df, polarity, filtered_ions, filtered_formulas);

  // ── neutral loss table ────────────────────────────────────────────────────
  Rcpp::NumericVector   nl_mass = neutral_loss_table["Accurate.Mass"];
  Rcpp::CharacterVector nl_name = neutral_loss_table.containsElementNamed("Formula")
    ? (Rcpp::CharacterVector)neutral_loss_table["Formula"]
  : (Rcpp::CharacterVector)neutral_loss_table["Neutral.Loss"];
  Rcpp::CharacterVector nl_pos  = neutral_loss_table["Pos"];
  Rcpp::CharacterVector nl_neg  = neutral_loss_table["Neg"];

  // ── iso label lookup map (from Phase 1 global run — not re-computed) ─────
  // Phase 2 features are not part of any library group, so there is no
  // in-group isotope run here.  All iso evidence in Phase 2 is pool-level
  // (_exp) by definition.  The "_exp" tag is appended when writing
  // iso_annotation so it is consistent with Phase 1 pool-fallback labels.
  std::unordered_map<std::string, std::string, KeyHash> iso_label_map;
  {
    Rcpp::CharacterVector il_fn    = iso_labels_df["feature_name"];
    Rcpp::CharacterVector il_label = iso_labels_df["iso_label"];
    for (int i = 0; i < il_fn.size(); ++i)
      iso_label_map[Rcpp::as<std::string>(il_fn[i])] =
        Rcpp::as<std::string>(il_label[i]);
  }

  // ── Phase 1 annotated set (skip these in pool scan) ───────────────────────
  std::unordered_set<std::string, KeyHash> p1_set;
  for (int i = 0; i < phase1_features.size(); ++i)
    p1_set.insert(Rcpp::as<std::string>(phase1_features[i]));

  // ── parse FT ─────────────────────────────────────────────────────────
  if (FT.nrows() == 0 || !FT.containsElementNamed("feature_name"))
    return Rcpp::List();

  Rcpp::CharacterVector pfn = FT["feature_name"];
  auto get_col = [&](const char* a, const char* b,
                     const char* c) -> Rcpp::NumericVector {
                       if (FT.containsElementNamed(a)) return FT[a];
                       if (FT.containsElementNamed(b)) return FT[b];
                       if (FT.containsElementNamed(c)) return FT[c];
                       return Rcpp::NumericVector(pfn.size(), NA_REAL);
                     };
  Rcpp::NumericVector pmz = get_col("mzmed", "mz",   "m/z");
  Rcpp::NumericVector prt = get_col("rtmed", "rt",   "retention_time");
  Rcpp::NumericVector pin = [&]() -> Rcpp::NumericVector {
    for (const char* nm : {"maxo","into","intb","intensity","intensity_median",
         "median_intensity","maxInt","intf"}) {
      if (FT.containsElementNamed(nm)) return FT[nm];
    }
    return Rcpp::NumericVector(pfn.size(), NA_REAL);
  }();
  int np = pfn.size();

  // ── parse confirmed_M_df ──────────────────────────────────────────────────
  Rcpp::NumericVector   cm_M       = confirmed_M_df["M"];
  Rcpp::CharacterVector cm_name    = confirmed_M_df["compound_name"];
  Rcpp::CharacterVector cm_lib     = confirmed_M_df["lib_group"];
  Rcpp::NumericVector   cm_ref_rt  = confirmed_M_df["ref_rt"];
  int n_M = cm_M.size();

  // ── Phase 2 expansion ─────────────────────────────────────────────────────
  std::vector<Rcpp::List> rows;

  for (int m = 0; m < n_M; ++m) {
    double M       = cm_M[m];
    double ref_rt  = cm_ref_rt[m];
    std::string compound_name = Rcpp::as<std::string>(cm_name[m]);
    std::string lib_group     = Rcpp::as<std::string>(cm_lib[m]);

    if (Rcpp::NumericVector::is_na(M)) continue;

    std::vector<double> theo_mz = parse_adduct_vector(M, filtered_formulas);
    size_t L = std::min(theo_mz.size(), filtered_ions.size());

    // Aggregate all matching adducts per pool feature into one row
    struct ExpFeatAgg {
      std::string fn;
      double mz, rt, intensity;
      std::vector<std::string> adduct_ions;
      std::vector<double>      ppm_errs, mda_errs;
    };
    std::unordered_map<std::string, ExpFeatAgg, KeyHash> exp_map;

    for (int p = 0; p < np; ++p) {
      std::string fn_p = trim(Rcpp::as<std::string>(pfn[p]));
      if (p1_set.count(fn_p)) continue;
      if (Rcpp::NumericVector::is_na(pmz[p])) continue;

      double feat_mz = pmz[p];
      double feat_rt = prt[p];

      if (!Rcpp::NumericVector::is_na(ref_rt) &&
          !Rcpp::NumericVector::is_na(feat_rt)) {
          if (std::fabs(feat_rt - ref_rt) > rt_exp_threshold) continue;
      }

      for (size_t k = 0; k < L; ++k) {
        double theo = theo_mz[k];
        if (std::isnan(theo)) continue;
        double tol      = mz_tol_mix(theo, mz_threshold, ppm_exp);
        double abs_diff = std::fabs(feat_mz - theo);
        if (abs_diff > tol) continue;
        double ppm_err = (theo > 0.0) ? (abs_diff / theo * 1e6) : NA_REAL;
        if (!std::isnan(ppm_err) && ppm_err > ppm_exp) continue;

        auto& agg = exp_map[fn_p];
        agg.fn        = fn_p;
        agg.mz        = feat_mz;
        agg.rt        = feat_rt;
        agg.intensity = Rcpp::NumericVector::is_na(pin[p]) ? NA_REAL : (double)pin[p];
        agg.adduct_ions.push_back(filtered_ions[k]);
        agg.ppm_errs.push_back(ppm_err);
        agg.mda_errs.push_back(abs_diff * 1e3);
      }
    }

    if (exp_map.empty()) continue;

    // Flatten → sort by mz for NL pass
    struct ExpFlat {
      std::string fn, adduct_str, ppm_str, mda_str, conf;
      double mz, rt, intensity, best_ppm;
    };
    std::vector<ExpFlat> exp_feats;
    exp_feats.reserve(exp_map.size());

    for (auto& kv : exp_map) {
      const ExpFeatAgg& agg = kv.second;
      std::vector<std::string> ann_p, ppm_p, mda_p;
      double best_ppm = std::numeric_limits<double>::quiet_NaN();
      double best_mda = std::numeric_limits<double>::quiet_NaN();
      for (size_t i = 0; i < agg.adduct_ions.size(); ++i) {
        std::ostringstream a, p, md;
        a << agg.adduct_ions[i] << " " << std::fixed << std::setprecision(6) << M;
        if (!std::isnan(agg.ppm_errs[i]))
          p << std::fixed << std::setprecision(2) << agg.ppm_errs[i] << " ppm";
        else p << "NA";
        md << std::fixed << std::setprecision(4) << agg.mda_errs[i] << " mDa";
        ann_p.push_back(a.str());
        ppm_p.push_back(p.str());
        mda_p.push_back(md.str());
        if (std::isnan(best_ppm) ||
            (!std::isnan(agg.ppm_errs[i]) && agg.ppm_errs[i] < best_ppm))
          best_ppm = agg.ppm_errs[i];
        if (std::isnan(best_mda) || agg.mda_errs[i] < best_mda)
          best_mda = agg.mda_errs[i];
      }
      std::string conf = std::isnan(best_ppm) ? "Level3_Tentative"
      : ppm_to_tier(best_ppm);
      exp_feats.push_back({
        agg.fn,
        join_annotations(ann_p),
        join_annotations(ppm_p),
        join_annotations(mda_p),
        conf,
        agg.mz, agg.rt, agg.intensity, best_ppm
      });
    }

    // Sort by mz — required by run_nl
    std::sort(exp_feats.begin(), exp_feats.end(),
              [](const ExpFlat& a, const ExpFlat& b){ return a.mz < b.mz; });

    int n_exp = (int)exp_feats.size();
    Rcpp::CharacterVector exp_fn_r(n_exp);
    Rcpp::NumericVector   exp_mz_r(n_exp), exp_rt_r(n_exp), exp_in_r(n_exp);
    for (int i = 0; i < n_exp; ++i) {
      exp_fn_r[i] = exp_feats[i].fn;
      exp_mz_r[i] = exp_feats[i].mz;
      exp_rt_r[i] = exp_feats[i].rt;
      exp_in_r[i] = exp_feats[i].intensity;
    }

    // Build iso label vector for Phase 2 features (mz-sorted, same order as
    // exp_mz_r) so run_nl can skip confirmed satellites as NL targets.
    std::vector<std::string> exp_iso_labels(n_exp, "");
    for (int i = 0; i < n_exp; ++i) {
      auto iit = iso_label_map.find(exp_feats[i].fn);
      if (iit != iso_label_map.end() && !iit->second.empty())
        exp_iso_labels[i] = iit->second;
    }

    // NL within expansion group (RT gate = rt_exp_threshold)
    std::vector<std::string> exp_nl = run_nl(
      exp_mz_r, exp_rt_r, exp_fn_r,
      nl_mass, nl_name, nl_pos, nl_neg,
      polarity, mz_threshold, nl_ppm, rt_exp_threshold,
      exp_iso_labels);

    // ── emit rows ─────────────────────────────────────────────────────────
    for (int i = 0; i < n_exp; ++i) {
      const ExpFlat& ef = exp_feats[i];

      // iso lookup from global map — always _exp in Phase 2 since there is
      // no library group to run a per-group iso against.
      std::string iso_str = "";
      {
        auto it2 = iso_label_map.find(ef.fn);
        if (it2 != iso_label_map.end() && !it2->second.empty()) {
          std::string raw = it2->second;
          size_t b1 = raw.find('['), e1 = raw.find(']');
          if (b1 != std::string::npos && e1 != std::string::npos) {
            std::string grp_num = lib_group;
            size_t us = grp_num.rfind('_');
            if (us != std::string::npos) grp_num = grp_num.substr(us + 1);
            raw = "[grp" + grp_num + "_" +
                  raw.substr(b1 + 1, e1 - b1 - 1) + "]" +
                  raw.substr(e1 + 1);
          }
          iso_str = raw + "_exp";
        }
      }

      // ── Confidence tier assignment for Phase 2 ────────────────────────
      // All Phase 2 features enter via adduct match → base tier is
      // Tier1a_exp_AnnoMatch (priority 6).
      //
      // NL evidence (exp_nl) is recorded in NL_match for traceability but
      // must NEVER overwrite reference_tier when an adduct match exists.
      // Tier1a_exp_AnnoMatch (priority 6) is stronger than
      // Tier1d_exp_NeutralLoss (priority 9) — NL only fires when there is
      // genuinely NO adduct match (ef.adduct_str empty).
      //
      // Pool-level isotope is recorded in iso_annotation but does NOT
      // change the tier — adduct match is already stronger evidence.
      //
      // Multimer order is parsed from the leading digit in the adduct ion
      // name (e.g. "3" in "[3M+2Na-H]+"); absence of a leading digit means
      // monomer (order 1).
      {
        double multimer_order = 1.0;
        if (!ef.adduct_str.empty()) {
          size_t bpos = ef.adduct_str.find('[');
          if (bpos != std::string::npos && bpos + 1 < ef.adduct_str.size() &&
              std::isdigit(static_cast<unsigned char>(ef.adduct_str[bpos + 1]))) {
            size_t dpos = bpos + 1;
            std::string digits;
            while (dpos < ef.adduct_str.size() &&
                   std::isdigit(static_cast<unsigned char>(ef.adduct_str[dpos]))) {
              digits += ef.adduct_str[dpos];
              ++dpos;
            }
            if (!digits.empty()) multimer_order = std::stod(digits);
          }
        }
        double effective_M = (!Rcpp::NumericVector::is_na(M))
          ? M * multimer_order
        : NA_REAL;
        if (!Rcpp::NumericVector::is_na(effective_M) &&
            ef.mz > effective_M * 1.5)
          iso_str = "";
      }

      // ── Suppress NL for Phase 2 isotope satellites ──────────────────────
      // Same logic as Phase 1: a C13 satellite cannot be a genuine NL
      // fragment. Clear exp_nl[i] if the global iso label marks this
      // feature as a satellite ([M+k] with k >= 1, any charge state).
      if (!iso_str.empty()) {
        // iso_str already has "_exp" suffix; strip for parsing
        std::string l = iso_str;
        size_t exp_pos = l.find("_exp");
        if (exp_pos != std::string::npos) l = l.substr(0, exp_pos);
        size_t e1 = l.find(']');
        size_t b2 = (e1 != std::string::npos) ? l.find('[', e1 + 1) : std::string::npos;
        if (b2 != std::string::npos) {
          std::string ion_part = l.substr(b2 + 1);
          if (ion_part.size() >= 2 && ion_part[0] == 'M' && ion_part[1] == '+')
            exp_nl[i] = "";
        }
      }

      std::string tier;
      if (!ef.adduct_str.empty()) {
        // Adduct match confirmed — always Tier1a_exp regardless of NL
        tier = "Tier1a_exp_AnnoMatch";
      } else if (!exp_nl[i].empty()) {
        // No adduct match but NL found — use NL tier
        tier = "Tier1d_exp_NeutralLoss";
      } else {
        // Should not happen: every exp_feat was selected by adduct scan,
        // but guard defensively.
        tier = "Tier1a_exp_AnnoMatch";
      }

      rows.push_back(Rcpp::List::create(
          Rcpp::Named("feature_name")          = ef.fn,
          Rcpp::Named("mz")                    = ef.mz,
          Rcpp::Named("rt")                    = ef.rt,
          Rcpp::Named("intensity")             = ef.intensity,
          Rcpp::Named("Compound_name")         = compound_name,
          Rcpp::Named("Lib_group")             = lib_group,
          Rcpp::Named("Compound_M")            = M,
          Rcpp::Named("PI_match")              = std::string(""),
          Rcpp::Named("grp_anno_match")      = std::string(""),
          Rcpp::Named("grp_anno_ppm")        = std::string(""),
          Rcpp::Named("grp_anno_mda")        = std::string(""),
          Rcpp::Named("grp_anno_confidence") = std::string(""),
          Rcpp::Named("exp_anno_match")      = ef.adduct_str,
          Rcpp::Named("exp_anno_ppm")        = ef.ppm_str,
          Rcpp::Named("exp_anno_mda")        = ef.mda_str,
          Rcpp::Named("exp_anno_confidence") = ef.conf,
          Rcpp::Named("NL_match")              = exp_nl[i],
                                                       Rcpp::Named("iso_annotation")        = iso_str,
                                                       Rcpp::Named("reference_tier")       = tier,
                                                       Rcpp::Named("annotation_phase")      = "Phase2",
                                                       Rcpp::Named("Match_peak_score")      = NA_REAL,
                                                       Rcpp::Named("inner_cosine")          = NA_REAL,
                                                       Rcpp::Named("outer_cosine")          = NA_REAL
      ));
    }

    // ── Isotope satellite expansion for this compound ──────────────────────
    // For every adduct-matched feature above that carries an iso label,
    // look up all other members of its iso group and emit them as additional
    // rows.  Isotope satellites are never found by adduct scan (they match no
    // adduct formula), so without this pass they are permanently missing.
    //
    // Strategy: build a reverse map {gidx → {fn, label}} from the iso_label_map
    // for all features seen in this expansion group.  For each adduct-matched
    // feature that has a non-empty iso label, look up the group and emit any
    // member not yet emitted for this compound.
    //
    // We use iso_label_map (the global run results passed from Phase 1) as the
    // source of truth for what belongs to each group.  mz/rt/intensity come
    // directly from FT (the same pool arrays already in scope: pmz/prt/pin).
    {
      // Build a quick fn→index map for pool lookup
      // (pfn/pmz/prt/pin are already in scope from the pool parse above)

      // Collect {gidx → members} for groups touched by this compound's
      // adduct-matched features.  Only groups that have at least one member
      // in exp_feats need expansion.
      std::unordered_map<std::string, std::vector<std::pair<std::string,std::string>>,
                         KeyHash> grp_members;  // gidx → {fn, label}

      // First pass: find which iso groups are referenced by this compound
      std::unordered_set<std::string, KeyHash> touched_gidx;
      for (int i = 0; i < n_exp; ++i) {
        auto it2 = iso_label_map.find(exp_feats[i].fn);
        if (it2 == iso_label_map.end() || it2->second.empty()) continue;
        const std::string& glbl = it2->second;
        size_t b1 = glbl.find('['), e1 = glbl.find(']');
        if (b1 == std::string::npos || e1 == std::string::npos) continue;
        touched_gidx.insert(glbl.substr(b1 + 1, e1 - b1 - 1));
      }

      if (!touched_gidx.empty()) {
        // Second pass: collect all iso_label_map members whose gidx is touched
        for (auto& kv2 : iso_label_map) {
          if (kv2.second.empty()) continue;
          size_t b1 = kv2.second.find('['), e1 = kv2.second.find(']');
          if (b1 == std::string::npos || e1 == std::string::npos) continue;
          std::string gidx2 = kv2.second.substr(b1 + 1, e1 - b1 - 1);
          if (!touched_gidx.count(gidx2)) continue;
          grp_members[gidx2].push_back({kv2.first, kv2.second});
        }

        // Build set of fn already emitted for this compound in this M-loop pass
        // (includes the adduct-matched rows just pushed above)
        std::unordered_set<std::string, KeyHash> already_emitted_this_M;
        for (int i = 0; i < n_exp; ++i)
          already_emitted_this_M.insert(exp_feats[i].fn);

        // Build a pool fn→index map for coordinate lookup
        // Iterating pfn each time is acceptable since n_exp is small
        std::unordered_map<std::string, int, KeyHash> pool_idx_map;
        for (int p2 = 0; p2 < np; ++p2)
          pool_idx_map[trim(Rcpp::as<std::string>(pfn[p2]))] = p2;

        // Third pass: emit missing satellites
        for (auto& gkv : grp_members) {
          for (auto& member : gkv.second) {
            const std::string& sat_fn  = member.first;
            const std::string& sat_lbl = member.second;

            if (already_emitted_this_M.count(sat_fn)) continue;  // already in output
            if (p1_set.count(sat_fn)) continue;  // annotated in Phase 1 — skip

            // RT gate: satellite must be within rt_exp_threshold of ref_rt
            auto pidx_it = pool_idx_map.find(sat_fn);
            if (pidx_it == pool_idx_map.end()) continue;  // not in pool
            int pidx = pidx_it->second;

            if (!Rcpp::NumericVector::is_na(ref_rt) &&
                !Rcpp::NumericVector::is_na(prt[pidx])) {
                if (std::fabs(prt[pidx] - ref_rt) > rt_exp_threshold) continue;
            }

            already_emitted_this_M.insert(sat_fn);

            rows.push_back(Rcpp::List::create(
                Rcpp::Named("feature_name")          = sat_fn,
                Rcpp::Named("mz")                    = (double)pmz[pidx],
                                                                  Rcpp::Named("rt")                    = (double)prt[pidx],
                                                                                                                    Rcpp::Named("intensity")             = Rcpp::NumericVector::is_na(pin[pidx])
                             ? NA_REAL : (double)pin[pidx],
                                                    Rcpp::Named("Compound_name")         = compound_name,
                                                    Rcpp::Named("Lib_group")             = lib_group,
                                                    Rcpp::Named("Compound_M")            = M,
                                                    Rcpp::Named("PI_match")              = std::string(""),
                                                    Rcpp::Named("grp_anno_match")      = std::string(""),
                                                    Rcpp::Named("grp_anno_ppm")        = std::string(""),
                                                    Rcpp::Named("grp_anno_mda")        = std::string(""),
                                                    Rcpp::Named("grp_anno_confidence") = std::string(""),
                                                    Rcpp::Named("exp_anno_match")      = std::string(""),
                                                    Rcpp::Named("exp_anno_ppm")        = std::string(""),
                                                    Rcpp::Named("exp_anno_mda")        = std::string(""),
                                                    Rcpp::Named("exp_anno_confidence") = std::string(""),
                                                    Rcpp::Named("NL_match")              = std::string(""),
                                                    Rcpp::Named("iso_annotation")        = [&]() {
                                                      std::string raw = sat_lbl;
                                                      size_t b1 = raw.find('['), e1 = raw.find(']');
                                                      if (b1 != std::string::npos && e1 != std::string::npos) {
                                                        std::string grp_num = lib_group;
                                                        size_t us = grp_num.rfind('_');
                                                        if (us != std::string::npos) grp_num = grp_num.substr(us + 1);
                                                        raw = "[grp" + grp_num + "_" +
                                                              raw.substr(b1 + 1, e1 - b1 - 1) + "]" +
                                                              raw.substr(e1 + 1);
                                                      }
                                                      return raw + "_exp";
                                                    }(),
                                                    Rcpp::Named("reference_tier")       = std::string("Tier1b_exp_Isotope"),
                                                    Rcpp::Named("annotation_phase")      = std::string("Phase2"),
                                                    Rcpp::Named("Match_peak_score")      = NA_REAL,
                                                    Rcpp::Named("inner_cosine")          = NA_REAL,
                                                    Rcpp::Named("outer_cosine")          = NA_REAL
            ));
          }
        }
      }
    }
  } // end M loop

  Rcpp::List out(rows.size());
  for (int i = 0; i < (int)rows.size(); ++i) out[i] = rows[i];
  return out;
}
