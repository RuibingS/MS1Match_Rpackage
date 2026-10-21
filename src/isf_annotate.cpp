#include <Rcpp.h>
#include <cmath>
#include <string>
#include <vector>
using namespace Rcpp;

// ============================================================
// isf_annotate_cpp
//
// For each parent feature (Tier-1d, has MS2 consensus):
//   For each candidate feature in pool_mz / pool_ids:
//     Guard: skip candidates whose mz > parent_mz + mz_threshold_isf
//       (fragment ions are physically sub-precursor).
//     Check if candidate mz matches ANY fragment ion in
//     the parent's fragment_mz vector within tolerance:
//       tol = max(mz_threshold, candidate_mz * ppm_isf * 1e-6)
//     If matched → record best (lowest ppm) hit.
//
// Called twice from R wrapper:
//   Phase 1 — pool_mz / pool_ids = intra-group co-members
//   Phase 2 — pool_mz / pool_ids = RT-filtered full pool
//             (Phase 1 features already excluded by R wrapper)
//
// Returns a DataFrame with one row per matched candidate.
// ============================================================

// [[Rcpp::export]]
DataFrame isf_annotate_cpp(
    NumericVector   pool_mz,          // candidate feature m/z values
    CharacterVector pool_ids,         // candidate feature_ids (same length)
    NumericVector   fragment_mz,      // parent MS2 fragment ion m/z values
    std::string     parent_id,        // parent feature_id (Tier-1d)
    double          parent_mz,        // parent precursor m/z
    double          parent_rt,        // parent RT in seconds (NA_real_ if absent)
    NumericVector   pool_rt,          // candidate RT values (NA allowed)
    double          ppm_isf           = 10.0,
    double          mz_threshold_isf  = 0.005
) {
  int n_pool  = pool_mz.size();
  int n_frags = fragment_mz.size();

  // --- output vectors ---
  std::vector<std::string> out_id;
  std::vector<double>      out_mz, out_rt;
  std::vector<double>      out_frag_mz, out_ppm, out_mda, out_rt_diff;

  out_id.reserve(n_pool);
  out_mz.reserve(n_pool);

  for (int i = 0; i < n_pool; ++i) {
    double fmz = pool_mz[i];

    // Fragment ions are always sub-precursor by definition.
    // Reject any candidate whose mz exceeds the parent precursor mz.
    // A small tolerance (mz_threshold_isf) is added to allow for cases
    // where the precursor mz in the MS2 spectrum differs slightly from
    // the MS1 apex mz of the parent feature.
    if (!NumericVector::is_na(parent_mz) &&
        fmz > parent_mz + mz_threshold_isf) continue;

    double tol = std::max(mz_threshold_isf, fmz * ppm_isf * 1e-6);

    double best_diff  = 1e9;
    int    best_j     = -1;

    for (int j = 0; j < n_frags; ++j) {
      double diff = std::fabs(fmz - fragment_mz[j]);
      if (diff <= tol && diff < best_diff) {
        best_diff = diff;
        best_j    = j;
      }
    }

    if (best_j < 0) continue;   // no match — skip

    double frt     = NumericVector::is_na(pool_rt[i]) ? NA_REAL : pool_rt[i];
    double rt_diff = (NumericVector::is_na(parent_rt) || NumericVector::is_na(frt))
                       ? NA_REAL
                       : frt - parent_rt;

    out_id.push_back(as<std::string>(pool_ids[i]));
    out_mz.push_back(fmz);
    out_rt.push_back(frt);
    out_frag_mz.push_back(fragment_mz[best_j]);
    out_ppm.push_back(best_diff / fmz * 1e6);
    out_mda.push_back(best_diff * 1000.0);
    out_rt_diff.push_back(rt_diff);
  }

  // --- build output DataFrame ---
  // parent columns are scalars — repeat for every matched row
  int n_out = (int)out_id.size();

  std::vector<std::string> out_parent_id(n_out, parent_id);
  std::vector<double>      out_parent_mz(n_out, parent_mz);
  std::vector<double>      out_parent_rt(n_out, parent_rt);

  return DataFrame::create(
    _["feature_name"]        = out_id,
    _["mz"]                  = out_mz,
    _["rt"]                  = out_rt,
    _["parent_feature"]      = out_parent_id,
    _["parent_mz"]           = out_parent_mz,
    _["parent_rt"]           = out_parent_rt,
    _["matched_fragment_mz"] = out_frag_mz,
    _["isf_ppm"]             = out_ppm,
    _["isf_mda"]             = out_mda,
    _["rt_diff"]             = out_rt_diff,
    _["stringsAsFactors"]    = false
  );
}
