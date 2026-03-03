// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

// ---------- helpers: ranking for Spearman ----------
static inline std::vector<double> rank_average(const std::vector<double>& x) {
  const size_t n = x.size();
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::stable_sort(idx.begin(), idx.end(),
                   [&](size_t a, size_t b){ return x[a] < x[b]; });

  std::vector<double> ranks(n);
  size_t i = 0;
  while (i < n) {
    size_t j = i + 1;
    while (j < n && x[idx[j]] == x[idx[i]]) ++j;
    double r = ((double)i + (double)j - 1.0) / 2.0 + 1.0; // avg rank, 1-based
    for (size_t k = i; k < j; ++k) ranks[idx[k]] = r;
    i = j;
  }
  return ranks;
}

// ---------- helpers: Pearson (C++ only) ----------
static inline double corr_pearson(const std::vector<double>& a,
                                  const std::vector<double>& b) {
  const size_t n = a.size();
  if (n != b.size() || n < 3) return NA_REAL;

  double mean_a = std::accumulate(a.begin(), a.end(), 0.0) / (double)n;
  double mean_b = std::accumulate(b.begin(), b.end(), 0.0) / (double)n;

  double num = 0.0, den_a = 0.0, den_b = 0.0;
  for (size_t i = 0; i < n; ++i) {
    double da = a[i] - mean_a;
    double db = b[i] - mean_b;
    num   += da * db;
    den_a += da * da;
    den_b += db * db;
  }
  if (den_a <= 0.0 || den_b <= 0.0) return NA_REAL;
  return num / std::sqrt(den_a * den_b);
}

// ---------- align + correlate (numeric only) ----------
static inline double compute_eic_correlation(
    const std::vector<double>& rt1,
    const std::vector<double>& int1,
    const std::vector<double>& rt2,
    const std::vector<double>& int2,
    double rt_EIC_thr,
    const std::string& method
) {
  std::vector<double> v1, v2;
  v1.reserve(rt1.size());
  v2.reserve(rt1.size());

  for (size_t i = 0; i < rt1.size(); ++i) {
    const double y1 = int1[i];
    if (std::isnan(y1) || y1 <= 0.0) continue;

    double best_diff = rt_EIC_thr + 1.0;
    int best_j = -1;
    for (size_t j = 0; j < rt2.size(); ++j) {
      const double y2 = int2[j];
      if (std::isnan(y2) || y2 <= 0.0) continue;
      double d = std::fabs(rt1[i] - rt2[j]);
      if (d <= rt_EIC_thr && d < best_diff) { best_diff = d; best_j = (int)j; }
    }
    if (best_j >= 0) {
      v1.push_back(std::log10(y1));
      v2.push_back(std::log10(int2[(size_t)best_j]));
    }
  }

  if (v1.size() < 3) return NA_REAL;

  if (method == "spearman") {
    std::vector<double> r1 = rank_average(v1);
    std::vector<double> r2 = rank_average(v2);
    return corr_pearson(r1, r2);
  }
  // default pearson
  return corr_pearson(v1, v2);
}

// ---------- main ----------
// [[Rcpp::export]]
DataFrame EIC_cor_mat_pairs(
     List feature_EIC_combined_result,
     DataFrame feature_pairs,
     double rt_EIC_thr = 1.0,
     std::string method = "pearson",
     bool parallel = true,
     int n_threads = 0
 ) {
   // Map feature name -> index in the input list
   CharacterVector names = feature_EIC_combined_result.names();

   struct Series { std::vector<double> rt, in; bool ok=false; };
   std::unordered_map<std::string, Series> cache;
   cache.reserve(names.size());

   for (int i = 0; i < names.size(); ++i) {
     std::string fname = as<std::string>(names[i]);
     List it = feature_EIC_combined_result[i];
     Series s;

     if (it.containsElementNamed("EIC_df")) {
       DataFrame df = as<DataFrame>(it["EIC_df"]);                // <- explicit cast
       if (df.containsElementNamed("rt") && df.containsElementNamed("intensity")) {
         s.rt = as< std::vector<double> >(df["rt"]);
         s.in = as< std::vector<double> >(df["intensity"]);
         s.ok = !(s.rt.empty() || s.in.empty() || s.rt.size() != s.in.size());
       }
     }
     cache.emplace(fname, std::move(s));
   }

   CharacterVector feat1 = feature_pairs["feature_1"];
   CharacterVector feat2 = feature_pairs["feature_2"];
   const int n_pairs = feat1.size();

   std::vector<std::string> out_f1(n_pairs), out_f2(n_pairs);
   std::vector<double> out_sim(n_pairs, NA_REAL);

#ifdef _OPENMP
   if (parallel && n_threads > 0) omp_set_num_threads(n_threads);
#endif

#ifdef _OPENMP
   if (parallel) {
#pragma omp parallel for schedule(dynamic)
     for (int k = 0; k < n_pairs; ++k) {
       std::string n1 = as<std::string>(feat1[k]);
       std::string n2 = as<std::string>(feat2[k]);
       double sim = NA_REAL;

       auto it1 = cache.find(n1);
       auto it2 = cache.find(n2);
       if (it1 != cache.end() && it2 != cache.end() && it1->second.ok && it2->second.ok) {
         sim = compute_eic_correlation(it1->second.rt, it1->second.in,
                                       it2->second.rt, it2->second.in,
                                       rt_EIC_thr, method);
       }
       out_sim[k] = sim; out_f1[k] = n1; out_f2[k] = n2;
     }
   } else
#endif
{
  for (int k = 0; k < n_pairs; ++k) {
    std::string n1 = as<std::string>(feat1[k]);
    std::string n2 = as<std::string>(feat2[k]);
    double sim = NA_REAL;

    auto it1 = cache.find(n1);
    auto it2 = cache.find(n2);
    if (it1 != cache.end() && it2 != cache.end() && it1->second.ok && it2->second.ok) {
      sim = compute_eic_correlation(it1->second.rt, it1->second.in,
                                    it2->second.rt, it2->second.in,
                                    rt_EIC_thr, method);
    }
    out_sim[k] = sim; out_f1[k] = n1; out_f2[k] = n2;
  }
}

   return DataFrame::create(
     _["feature_1"]      = out_f1,
     _["feature_2"]      = out_f2,
     _["EIC_similarity"] = out_sim
   );
 }
