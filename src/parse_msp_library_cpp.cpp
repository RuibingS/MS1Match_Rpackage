#include <Rcpp.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <cctype>
#include <cmath>

using namespace Rcpp;

// =============================================================================
// Monoisotopic mass table
// =============================================================================
static const std::unordered_map<std::string, double> MONO_MASS = {
  {"H",   1.00782503207}, {"C",  12.00000000000}, {"N",  14.00307400480},
  {"O",  15.99491461956}, {"S",  31.97207100000}, {"P",  30.97376163000},
  {"F",  18.99840322000}, {"Cl", 34.96885268000}, {"Br", 78.91833710000},
  {"I", 126.90447300000}, {"Si", 27.97692653000}, {"B",  11.00930540000},
  {"Se", 79.91652100000}, {"Fe", 55.93493750000}, {"Cu", 62.92959750000},
  {"Zn", 63.92914220000}, {"Ca", 39.96259098000}, {"Mg", 23.98504170000},
  {"Na", 22.98976928000}, {"K",  38.96370668000}, {"Li",  6.01512279000},
  {"Co", 58.93319500000}, {"Mn", 54.93804500000}, {"Mo", 97.90540700000},
  {"Ni", 57.93534600000}, {"As", 74.92159700000}, {"Hg",201.97062600000},
  {"Cd",113.90335800000}, {"Pb",203.97304100000}, {"Al", 26.98153860000},
  {"Au",196.96656900000}, {"Ag",106.90509700000}
};

static double formula_to_exact_mass(const std::string &formula) {
  double mass = 0.0;
  size_t i = 0, n = formula.size();
  while (i < n) {
    if (!std::isupper((unsigned char)formula[i])) return -1.0;
    std::string elem;
    elem += formula[i++];
    while (i < n && std::islower((unsigned char)formula[i])) elem += formula[i++];
    int count = 0;
    while (i < n && std::isdigit((unsigned char)formula[i]))
      count = count * 10 + (formula[i++] - '0');
    if (count == 0) count = 1;
    auto it = MONO_MASS.find(elem);
    if (it == MONO_MASS.end()) return -1.0;
    mass += it->second * count;
  }
  return mass;
}

// [[Rcpp::export]]
Rcpp::NumericVector calculate_exact_mass_cpp(Rcpp::CharacterVector formulas) {
  int n = formulas.size();
  Rcpp::NumericVector out(n, NA_REAL);
  for (int i = 0; i < n; ++i) {
    if (formulas[i] == NA_STRING) continue;
    double m = formula_to_exact_mass(Rcpp::as<std::string>(formulas[i]));
    if (m > 0.0) out[i] = m;
  }
  return out;
}

// =============================================================================
// String helpers
// =============================================================================
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                  [](unsigned char c){ return !std::isspace(c); }));
}
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char c){ return !std::isspace(c); }).base(), s.end());
}
static std::string to_lower(const std::string &s) {
  std::string o = s;
  std::transform(o.begin(), o.end(), o.begin(), ::tolower);
  return o;
}
static std::string normalize_ion_mode(const std::string &s) {
  std::string v = to_lower(s);
  if (v=="positive"||v=="pos"||v=="p") return "POSITIVE";
  if (v=="negative"||v=="neg"||v=="n") return "NEGATIVE";
  return s;
}
// "msms" accepted as an alias for "MS2" in the
// e.g. "Spectrum_type: MSMS" to "MS2" like
// "Spectrum_type: MS2" 
static std::string normalize_spectrum_type(const std::string &s) {
  std::string v = to_lower(s);
  if (v=="ms"  || v=="ms1")  return "MS1";
  if (v=="ms2" || v=="msms") return "MS2";
  return s;
}


// =============================================================================
struct CompoundRecord {
  std::string name, formula, ion_mode, spectrum_type;
  std::string instrument_type, instrument, ionization_method;
  std::string precursor_mz, precursor_type, collision_energy;
  std::string retention_time, comments;
  std::string inchikey, inchi, smiles, db_num, splash;
  double exact_mass = NA_REAL;
  int    mw         = NA_INTEGER;
  std::vector<std::pair<std::string,std::string>> extra_fields;
  std::vector<double> mzs;
  std::vector<double> intensities;
};

// =============================================================================
// R list
// =============================================================================
static Rcpp::List record_to_rlist(const CompoundRecord &r) {
  std::vector<std::string> fnames;
  std::vector<SEXP>        fvalues;
  int n_fixed = 19;
  int n_extra = (int)r.extra_fields.size();
  fnames.reserve(n_fixed + n_extra + 1);
  fvalues.reserve(n_fixed + n_extra + 1);

  auto add_s = [&](const char* nm, const std::string &v) {
    if (!v.empty()) { fnames.push_back(nm); fvalues.push_back(Rcpp::wrap(v)); }
  };
  auto add_d = [&](const char* nm, double v) {
    if (!std::isnan(v)) { fnames.push_back(nm); fvalues.push_back(Rcpp::wrap(v)); }
  };
  auto add_i = [&](const char* nm, int v) {
    if (v != NA_INTEGER) { fnames.push_back(nm); fvalues.push_back(Rcpp::wrap(v)); }
  };

  add_s("Name",              r.name);
  add_s("Formula",           r.formula);
  add_d("ExactMass",         r.exact_mass);
  add_i("MW",                r.mw);
  add_s("Ion_mode",          r.ion_mode);
  add_s("Spectrum_type",     r.spectrum_type);
  add_s("Instrument_type",   r.instrument_type);
  add_s("Instrument",        r.instrument);
  add_s("Ionization_method", r.ionization_method);
  add_s("PrecursorMZ",       r.precursor_mz);
  add_s("Precursor_type",    r.precursor_type);
  add_s("Collision_energy",  r.collision_energy);
  add_s("RetentionTime",     r.retention_time);
  add_s("InChIKey",          r.inchikey);
  add_s("InChI",             r.inchi);
  add_s("SMILES",            r.smiles);
  add_s("DB#",               r.db_num);
  add_s("Splash",            r.splash);
  add_s("Comments",          r.comments);

  for (auto &kv : r.extra_fields) {
    fnames.push_back(kv.first);
    fvalues.push_back(Rcpp::wrap(kv.second));
  }

  // Peaks DataFrame
  fnames.push_back("Peaks");
  fvalues.push_back(Rcpp::DataFrame::create(
      Rcpp::Named("mz")            = Rcpp::NumericVector(r.mzs.begin(),
                  r.mzs.end()),
                  Rcpp::Named("rel_intensity") = Rcpp::NumericVector(r.intensities.begin(),
                              r.intensities.end())
  ));

  int nf = (int)fnames.size();
  Rcpp::List compound(nf);
  Rcpp::CharacterVector nms(nf);
  for (int k = 0; k < nf; ++k) {
    compound[k] = fvalues[k];
    nms[k]      = fnames[k];
  }
  compound.attr("names") = nms;
  return compound;
}

// =============================================================================
// Internal parser
// =============================================================================
static Rcpp::List parse_msp_internal(
    std::string file_path,
    Rcpp::Nullable<std::string> spectrum_type,
    Rcpp::Nullable<std::string> ionPolarity,
    Rcpp::Nullable<std::string> ioniMethod,
    Rcpp::Nullable<std::string> colEnergy,
    Rcpp::Nullable<std::string> instName,
    std::string dedup_by,
    int reserve_n,
    int max_records) {

  std::ifstream infile(file_path, std::ios::in);
  if (!infile.is_open()) Rcpp::stop("Could not open file: " + file_path);

  bool f_spec = spectrum_type.isNotNull();
  bool f_pol  = ionPolarity.isNotNull();
  bool f_meth = ioniMethod.isNotNull();
  bool f_en   = colEnergy.isNotNull();
  bool f_inst = instName.isNotNull();
  std::string s_spec = f_spec ? normalize_spectrum_type(
    Rcpp::as<std::string>(spectrum_type)) : "";
  std::string s_pol  = f_pol  ? normalize_ion_mode(
    Rcpp::as<std::string>(ionPolarity))   : "";
  std::string s_meth = f_meth ? Rcpp::as<std::string>(ioniMethod) : "";
  std::string s_en   = f_en   ? Rcpp::as<std::string>(colEnergy)  : "";
  std::string s_inst = f_inst ? Rcpp::as<std::string>(instName)   : "";

  std::vector<CompoundRecord> all_records;
  all_records.reserve(reserve_n);

  CompoundRecord cur;
  cur.mzs.reserve(128);
  cur.intensities.reserve(128);

  bool in_peaks = false, has_msms = false;
  bool valid = true, has_em = false, has_any = false;
  int  n_exp = 0, n_read = 0;
  int  total_count = 0, filtered_count = 0;

  // ---- flush: finalize current compound ----
  auto flush = [&]() {
    if (!has_any || n_read == 0) { if (has_any) filtered_count++; return; }
    total_count++;

    // infer spectrum_type before the filter check uses it

    if (cur.spectrum_type.empty())
      cur.spectrum_type = has_msms ? "MS2" : "MS1";

    bool accept = valid;
    if (f_spec && accept) {
      bool is_ms2 = has_msms || cur.spectrum_type == "MS2";
      if (s_spec == "MS2" && !is_ms2) accept = false;
      if (s_spec == "MS1" &&  is_ms2) accept = false;
    }
    if (!accept) { filtered_count++; return; }

    // Normalize to base peak = 100
    double mx = 0.0;
    for (double v : cur.intensities) if (v > mx) mx = v;
    if (mx > 0.0)
      for (double &v : cur.intensities) v = v / mx * 100.0;

    all_records.push_back(std::move(cur));
  };

  auto reset = [&]() {
    cur = CompoundRecord();
    cur.mzs.reserve(128);
    cur.intensities.reserve(128);
    in_peaks = has_msms = has_em = has_any = false;
    valid = true;
    n_exp = n_read = 0;
  };

  std::string line;
  line.reserve(2048);

  while (std::getline(infile, line)) {

    if (max_records > 0 && (int)all_records.size() >= max_records) break;

    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line.front() == ' ') ltrim(line);
    if (!line.empty() && line.back()  == ' ') rtrim(line);

    // Blank line = end of record
    if (line.empty()) {
      if (has_any) { flush(); reset(); }
      continue;
    }

    // New Name: = start of new record
    if (has_any &&
        line.size() > 5 &&
        line[0]=='N' && line[1]=='a' && line[2]=='m' &&
        line[3]=='e' && line[4]==':') {
      flush(); reset();
    }

    // Peak data lines
    if (in_peaks && n_read < n_exp) {
      char c0 = line[0];
      if (std::isdigit((unsigned char)c0) || c0 == '.') {
        std::istringstream iss(line);
        double mz, inten;
        while (iss >> mz >> inten && n_read < n_exp) {
          cur.mzs.push_back(mz);
          cur.intensities.push_back(inten);
          ++n_read;
        }
        continue;
      } else if (c0 != '-') {
        in_peaks = false;
      }
    }

    // Key : Value
    size_t colon = line.find(':');
    if (colon == std::string::npos) continue;

    size_t klen = colon;
    while (klen > 0 && line[klen-1] == ' ') --klen;
    std::string key(line.c_str(), klen);

    size_t vs = colon + 1;
    while (vs < line.size() && line[vs] == ' ') ++vs;
    std::string value = line.substr(vs);
    while (!value.empty() && std::isspace((unsigned char)value.back()))
      value.pop_back();

    std::string kl = to_lower(key);
    has_any = true;

    if      (kl == "name")             { cur.name = value; }
    else if (kl == "formula")          {
      cur.formula = value;
      if (!has_em) {
        double em = formula_to_exact_mass(value);
        if (em > 0.0) cur.exact_mass = em;
      }
    }
    else if (kl == "exactmass")        {
      try { cur.exact_mass = std::stod(value); has_em = true; } catch(...) {}
    }
    else if (kl == "mw")               {
      try { cur.mw = std::stoi(value); } catch(...) {}
    }
    else if (kl == "ion_mode" || kl == "ionpolarity") {
      cur.ion_mode = normalize_ion_mode(value);
      if (f_pol && cur.ion_mode != s_pol) valid = false;
    }
    else if (kl == "spectrum_type")    { cur.spectrum_type = normalize_spectrum_type(value); }
   
    else if (kl == "msms" || kl == "ms2") { has_msms = true; }
    else if (kl == "precursor_type")   { cur.precursor_type = value; }
    else if (kl == "precursormz" || kl == "preion") { cur.precursor_mz = value; }
    else if (kl == "instrument_type")  {
      cur.instrument_type = value;
      if (f_meth && value != s_meth) valid = false;
    }
    else if (kl == "instrument")       {
      cur.instrument = value;
      if (f_inst && value != s_inst) valid = false;
    }
    else if (kl == "ionization_method" || kl == "ionimethod") {
      cur.ionization_method = value;
    }
    else if (kl == "collision_energy") {
      cur.collision_energy = value;
      if (f_en && value != s_en) valid = false;
    }
    else if (kl == "rettime" || kl == "retentiontime") { cur.retention_time = value; }
    else if (kl == "comments" || kl == "comment")      { cur.comments = value; }
    else if (kl == "inchikey")  { cur.inchikey = value; }
    else if (kl == "inchi")     { cur.inchi    = value; }
    else if (kl == "smiles")    { cur.smiles   = value; }
    else if (kl == "db#")       { cur.db_num   = value; }
    else if (kl == "splash")    { cur.splash   = value; }
    else if (kl == "num peaks") {
      try { n_exp = std::stoi(value); } catch(...) { n_exp = 0; }
      in_peaks = (n_exp > 0);
    }
    else { cur.extra_fields.emplace_back(std::move(key), std::move(value)); }
  }

  if (has_any) { flush(); reset(); }
  infile.close();

  Rcpp::Rcout << "[parser] Parsed:          " << total_count        << "\n"
              << "[parser] Filtered:        " << filtered_count     << "\n"
              << "[parser] Pre-dedup total: " << all_records.size() << "\n";

  // =============================================================================
  // Deduplication
  // =============================================================================
  std::vector<int> keep_idx;

  if (dedup_by == "none") {
    keep_idx.resize(all_records.size());
    std::iota(keep_idx.begin(), keep_idx.end(), 0);
  } else {
    auto make_key = [&](const CompoundRecord &r) -> std::string {
      if (dedup_by == "name_only")        return r.name;
      if (dedup_by == "name_adduct_mode") return r.name+"|"+r.precursor_type+"|"+r.ion_mode;
      return r.name+"|"+r.precursor_type+"|"+r.ion_mode+"|"+r.collision_energy;
    };
    std::unordered_map<std::string,int> best;
    best.reserve(all_records.size());
    for (int i = 0; i < (int)all_records.size(); ++i) {
      std::string k = make_key(all_records[i]);
      auto it = best.find(k);
      if (it == best.end()) {
        best[k] = i;
      } else if (all_records[i].mzs.size() > all_records[it->second].mzs.size()) {
        it->second = i;
      }
    }
    keep_idx.reserve(best.size());
    for (auto &kv : best) keep_idx.push_back(kv.second);
  }

  int n_out     = (int)keep_idx.size();
  int dup_count = (int)all_records.size() - n_out;

  Rcpp::Rcout << "[parser] Duplicates:      " << dup_count << "\n"
              << "[parser] Building result: " << n_out << " compounds...\n";

  Rcpp::List result(n_out);
  for (int ri = 0; ri < n_out; ++ri)
    result[ri] = record_to_rlist(all_records[keep_idx[ri]]);

  Rcpp::Rcout << "[parser] Done.\n";
  return result;
}


// [[Rcpp::export]]
Rcpp::List parse_msp_library_cpp(
    std::string file_path,
    Rcpp::Nullable<std::string> spec = R_NilValue,
    Rcpp::Nullable<std::string> pol  = R_NilValue,
    Rcpp::Nullable<std::string> meth = R_NilValue,
    Rcpp::Nullable<std::string> cen  = R_NilValue,
    Rcpp::Nullable<std::string> inst = R_NilValue,
    std::string dedup  = "name_adduct_mode_energy",
    int reserven       = 60000,
    int maxrec         = -1) {
  return parse_msp_internal(file_path, spec, pol, meth, cen, inst,
                            dedup, reserven, maxrec);
}
