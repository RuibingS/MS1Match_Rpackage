#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <limits>
#include "utils_check.h"
#include "utils_sort.h"
using namespace Rcpp;


// --------------------------------------------------------
// PImatch_target_rcpp
//
// Matches feature table m/z values against a normalised
// target list (output of normalize_target_list()).
//
// Standard input columns expected from Comp_data:
//   Name       — compound name
//   exactMass  — monoisotopic / neutral exact mass
//   [M+H]+     — positive primary ion m/z
//   [M-H]-     — negative primary ion m/z
//   RT         — retention time (numeric, NA if absent)
//
// Arguments
//   FT           feature table  (needs feature_name, mzmed/mz)
//   Comp_data    normalised target list (see above)
//   polarity     "pos" or "neg"
//   diff_mz_thr  absolute m/z tolerance (Da)
//   ppm          ppm m/z tolerance (applied as max of ppm and diff_mz_thr)
//   rt_thr_compound  RT window (seconds) for feature vs compound RT;
//                    matching is skipped when compound RT is NA
// --------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List PImatch_target_rcpp(Rcpp::DataFrame FT,
                               Rcpp::DataFrame Comp_data,
                               std::string     polarity,
                               double          diff_mz_thr,
                               int             ppm,
                               double          rt_thr_compound) {

  int nrow_FT   = FT.nrow();
  int nrow_Comp = Comp_data.nrow();

  // ---- feature table columns ----
  NumericVector   FT_mz  = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  CharacterVector FT_name = get_required_column<CharacterVector>(FT, "feature_name");

  NumericVector FT_rt;
  bool has_ft_rt = false;
  try {
    FT_rt      = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention_time");
    has_ft_rt  = (FT_rt.size() == FT_mz.size());
  } catch (...) {
    has_ft_rt = false;
  }

  // ---- compound table columns (normalised schema) ----
  CharacterVector Comp_name     = Comp_data["Name"];
  NumericVector   Comp_mass     = Comp_data["exactMass"];

  // Choose ion m/z column based on polarity
  NumericVector comp_mz;
  std::string   ion_label;
  if (polarity == "pos") {
    comp_mz   = Comp_data["[M+H]+"];
    ion_label = "[M+H]+";
  } else if (polarity == "neg") {
    comp_mz   = Comp_data["[M-H]-"];
    ion_label = "[M-H]-";
  } else {
    stop("Invalid polarity: must be 'pos' or 'neg'");
  }

  // Optional compound RT column
  NumericVector Comp_rt;
  bool has_comp_rt = false;
  try {
    Comp_rt      = Comp_data["RT"];
    has_comp_rt  = true;
  } catch (...) {
    has_comp_rt = false;
  }

  // ---- matching ----
  Rcpp::List output(nrow_FT);

  for (int i = 0; i < nrow_FT; ++i) {

    double fmz = FT_mz[i];
    double frt = has_ft_rt ? FT_rt[i] : NA_REAL;

    CharacterVector feat_names_out;
    NumericVector   feat_mz_out;
    NumericVector   feat_rt_out;
    CharacterVector comp_names_out;
    CharacterVector ion_label_out;
    NumericVector   exact_mass_out;
    NumericVector   comp_mz_out;
    NumericVector   comp_rt_out;
    NumericVector   mz_diff_out;
    NumericVector   ppm_err_out;

    for (int j = 0; j < nrow_Comp; ++j) {

      double cmz = comp_mz[j];

      // Skip compound rows where the ion m/z or exact mass is NA/NaN.

      if (ISNAN(cmz) || ISNAN(Comp_mass[j])) continue;

      double mz_tol  = std::max(diff_mz_thr, cmz * ppm * 1e-6);
      double mz_diff = std::fabs(fmz - cmz);

      if (mz_diff > mz_tol) continue;

      // RT gate: only applied when:
      //   1. rt_thr_compound is a valid (non-NA/NaN) threshold — if the user
      //      passes NA or NULL from R it arrives as NaN; NaN comparisons are
      //      always false in C++ so the gate would silently never reject.
      //   2. Both the feature RT and the compound RT are non-NA.
      double crt = (has_comp_rt && !ISNA(Comp_rt[j])) ? Comp_rt[j] : NA_REAL;
      if (!ISNAN(rt_thr_compound) &&
          has_ft_rt && !ISNA(frt) && !ISNA(crt)) {
        if (std::fabs(frt - crt) > rt_thr_compound) continue;
      }

      feat_names_out.push_back(as<std::string>(FT_name[i]));
      feat_mz_out.push_back(fmz);
      feat_rt_out.push_back(has_ft_rt ? frt : NA_REAL);
      comp_names_out.push_back(as<std::string>(Comp_name[j]));
      ion_label_out.push_back(ion_label + " " + std::to_string(Comp_mass[j]));
      exact_mass_out.push_back(Comp_mass[j]);
      comp_mz_out.push_back(cmz);
      comp_rt_out.push_back(ISNA(crt) ? NA_REAL : crt);
      mz_diff_out.push_back(mz_diff);
      ppm_err_out.push_back(1e6 * mz_diff / cmz);
    }

    output[i] = List::create(
      _["Feature_name"] = feat_names_out,
      _["Feature_mz"]   = feat_mz_out,
      _["Feature_rt"]   = feat_rt_out,
      _["Comp_name"]    = comp_names_out,
      _["PI_name"]      = ion_label_out,
      _["exactMass"]    = exact_mass_out,
      _["Comp_mz"]      = comp_mz_out,
      _["Comp_rt"]      = comp_rt_out,
      _["mz_diff"]      = mz_diff_out,
      _["ppm_error"]    = ppm_err_out
    );
  }

  return output;
}
