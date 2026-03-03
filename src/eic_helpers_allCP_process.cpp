// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <cstdio>     // std::snprintf
#include <cmath>

using namespace Rcpp;

#ifdef _OPENMP
#include <omp.h>
#endif

// ---------- helpers: detect matrix-like vs list-of-rows --------------------
static inline bool is_matrix_like(const List& x) {
  return x.hasAttribute("dim");
}
static inline S4 get_cell_matrix_like(const List& x, int i, int j, int nrows) {
  // column-major index
  return as<S4>( x[j * nrows + i] );
}
static inline S4 get_cell_rowlist(const List& x, int i, int j) {
  List row = x[i];
  return as<S4>( row[j] );
}
static inline S4 get_cell_flex(const List& x, int i, int j, int nrows) {
  return is_matrix_like(x) ? get_cell_matrix_like(x, i, j, nrows)
    : get_cell_rowlist(x, i, j);
}

// ---------- helpers: chromPeaks column access and row -> data.frame --------
static inline int col_index(const NumericMatrix& m, const std::string& name) {
  List dimn = m.attr("dimnames");
  if (dimn.size() != 2) {
    return -1;
  }
  CharacterVector cn = dimn[1];
  for (int k = 0; k < cn.size(); ++k) {
    if (name == as<std::string>(cn[k])) {
      return k;
    }
  }
  return -1;
}

static inline DataFrame cp_row_df(const NumericMatrix& pks, int r) {
  auto get = [&](const std::string& nm) {
    int c = col_index(pks, nm);
    return (c >= 0) ? NumericVector::create(pks(r, c))
      : NumericVector::create(NA_REAL);
  };
  return DataFrame::create(
    Named("mz")     = get("mz"),
    Named("mzmin")  = get("mzmin"),
    Named("mzmax")  = get("mzmax"),
    Named("rt")     = get("rt"),
    Named("rtmin")  = get("rtmin"),
    Named("rtmax")  = get("rtmax"),
    Named("into")   = get("into"),
    Named("intb")   = get("intb"),
    Named("maxo")   = get("maxo"),
    Named("sn")     = get("sn"),
    Named("sample") = (col_index(pks,"sample")>=0 ? get("sample") :
                         (col_index(pks,"column")>=0 ? get("column") :
                            NumericVector::create(NA_REAL)))
  );
}

// ---------- subset EIC by [rtmin, rtmax] (+ pad); drop NA intensities -------
static inline DataFrame subset_eic(const NumericVector& rt,
                                   const NumericVector& inten,
                                   double lo, double hi, double pad) {
  double a = lo - pad, b = hi + pad;
  std::vector<double> out_rt, out_in;
  out_rt.reserve(rt.size());
  out_in.reserve(inten.size());
  for (int k = 0; k < rt.size(); ++k) {
    double t = rt[k];
    if (NumericVector::is_na(t)) continue;
    if (t >= a && t <= b) {
      double y = inten[k];
      if (!NumericVector::is_na(y)) {
        out_rt.push_back(t);
        out_in.push_back(y);
      }
    }
  }
  return DataFrame::create(Named("rt") = out_rt, Named("intensity") = out_in);
}

// ---------- pure C++ snapshot of chromPeaks columns for safe parallel use ---
struct PksCols {
  int n = 0;
  bool has_rt=false, has_into=false, has_maxo=false, has_sn=false, has_rtmin=false, has_rtmax=false;
  std::vector<double> rt, into, maxo, sn, rtmin, rtmax;
};

static inline PksCols copy_pks_cols(const NumericMatrix& pks) {
  PksCols P;
  int n = pks.nrow(); P.n = n;
  if (n == 0) return P;

  auto take = [&](const char* nm, std::vector<double>& dst, bool& has){
    int c = col_index(pks, nm);
    if (c >= 0) {
      has = true; dst.resize(n);
      for (int r = 0; r < n; ++r) dst[r] = pks(r, c);
    }
  };

  take("rt",    P.rt,    P.has_rt);
  take("into",  P.into,  P.has_into);
  take("maxo",  P.maxo,  P.has_maxo);
  take("sn",    P.sn,    P.has_sn);
  take("rtmin", P.rtmin, P.has_rtmin);
  take("rtmax", P.rtmax, P.has_rtmax);

  return P;
}

// ---------- choose 1 CP row using snapshot columns (thread-safe) ------------
static inline int choose_cp_row_cols(const PksCols& P,
                                     const std::string& policy,
                                     double rtmed_i,
                                     double sn_min) {
  if (P.n == 0) return -1;

  // S/N prefilter
  std::vector<int> cand; cand.reserve(P.n);
  for (int r = 0; r < P.n; ++r) {
    if (P.has_sn && sn_min > 0 && !Rcpp::NumericVector::is_na(P.sn[r])) {
      if (P.sn[r] < sn_min) continue;
    }
    cand.push_back(r);
  }
  if (cand.empty()) return -1;

  auto pick_max = [&](const std::vector<double>& v, bool has)->int{
    if (!has || v.empty()) return -1;
    int best = cand[0]; double bestv = v[best];
    for (size_t k = 1; k < cand.size(); ++k) {
      int r = cand[k]; double z = v[r];
      if (Rcpp::NumericVector::is_na(z)) continue;
      if (Rcpp::NumericVector::is_na(bestv) || z > bestv) { bestv = z; best = r; }
    }
    return best;
  };

  if (policy == "max_into")  return pick_max(P.into, P.has_into);
  if (policy == "max_maxo")  return pick_max(P.maxo, P.has_maxo);

  // "medret" (default). Fallback to into -> maxo if missing.
  if (!P.has_rt || Rcpp::NumericVector::is_na(rtmed_i)) {
    int r = pick_max(P.into, P.has_into);
    if (r >= 0) return r;
    return pick_max(P.maxo, P.has_maxo);
  }
  int best = cand[0];
  double bestd = std::abs(P.rt[best] - rtmed_i);
  for (size_t k = 1; k < cand.size(); ++k) {
    int r = cand[k];
    double d = std::abs(P.rt[r] - rtmed_i);
    if (d < bestd) { bestd = d; best = r; }
  }
  return best;
}

// ---------- main ------------------------------------------------------------


// [[Rcpp::export]]
 List feature_EIC_process_full_fun_cpp(List input_list, int nrows, int ncols,
                                  CharacterVector eic_row_names,
                                  DataFrame ft,
                                  std::string policy = "medret",
                                  double sn_min = 0.0,
                                  double pad = 0.0,
                                  bool parallel = false,
                                  int n_threads = 0) {

   // check ft columns
   if (!ft.containsElementNamed("feature_name") || !ft.containsElementNamed("rtmed")) {
     stop("ft must contain columns 'feature_name' and 'rtmed'.");
   }
   CharacterVector ft_names = ft["feature_name"];
   NumericVector   ft_rtmed = ft["rtmed"];

   // build lookup map: feature_name -> rtmed
   std::unordered_map<std::string, double> rtmed_map;
   rtmed_map.reserve(ft_names.size());
   for (int i = 0; i < ft_names.size(); ++i) {
     std::string key = as<std::string>(ft_names[i]);
     if (rtmed_map.find(key) == rtmed_map.end()) {
       rtmed_map.emplace(key, NumericVector::is_na(ft_rtmed[i]) ? NA_REAL
                           : (double)ft_rtmed[i]);
     }
   }

   // feature names: provided or auto "FT%04d"
   CharacterVector feat_names(nrows);
   if (eic_row_names.size() == nrows) {
     feat_names = clone(eic_row_names);
   } else {
     for (int i = 0; i < nrows; ++i) {
       char buf[16];
       std::snprintf(buf, sizeof(buf), "FT%04d", i + 1);
       feat_names[i] = std::string(buf);
     }
   }

   // sample list names
   CharacterVector nm_eic(ncols), nm_cp(ncols);
   for (int j = 0; j < ncols; ++j) {
     nm_eic[j] = "EIC_df_sample_" + std::to_string(j + 1);
     nm_cp[j]  = "CP_sample_"     + std::to_string(j + 1);
   }

   const bool matrix_like = is_matrix_like(input_list);
   List out(nrows);

   // optional: set OMP threads
#ifdef _OPENMP
   if (parallel && n_threads > 0) omp_set_num_threads(n_threads);
#endif

   for (int i = 0; i < nrows; ++i) {
     std::string fname = as<std::string>(feat_names[i]);
     double rtmed_i = NA_REAL;
     auto it = rtmed_map.find(fname);
     if (it != rtmed_map.end()) {
       rtmed_i = it->second;
     }

     // collect per-sample CP matrices for chrom_peaks_df (sequential, R API)
     List chrom_peaks_list(ncols);
     for (int j = 0; j < ncols; ++j) {
       S4 xChrom = get_cell_flex(input_list, i, j, nrows);
       NumericMatrix pks = as<NumericMatrix>(xChrom.slot("chromPeaks"));
       chrom_peaks_list[j] = pks;
     }
     DataFrame chrom_peaks_df = Rcpp::as<DataFrame>(
       Rcpp::wrap(Rcpp::Language("do.call", "rbind", chrom_peaks_list).eval())
     );

     // snapshot pks columns for safe parallel selection
     std::vector<PksCols> Pks(ncols);
     for (int j = 0; j < ncols; ++j) {
       S4 xChrom = get_cell_flex(input_list, i, j, nrows);
       NumericMatrix pks = as<NumericMatrix>(xChrom.slot("chromPeaks"));
       Pks[j] = copy_pks_cols(pks);
     }

     // parallel-safe selection results
     std::vector<int>    chosen_idx(ncols, -1);
     std::vector<double> lo(ncols, NA_REAL), hi(ncols, NA_REAL);

     // -------- parallel (numeric-only) region --------------------------------
#ifdef _OPENMP
     if (parallel) {
#pragma omp parallel for schedule(dynamic)
       for (int j = 0; j < ncols; ++j) {
         const PksCols& P = Pks[j];
         int idx = choose_cp_row_cols(P, policy, rtmed_i, sn_min);
         chosen_idx[j] = idx;
         if (idx >= 0 && idx < P.n && P.has_rtmin && P.has_rtmax) {
           lo[j] = P.rtmin[idx];
           hi[j] = P.rtmax[idx];
         } else {
           chosen_idx[j] = -1;
         }
       }
     } else
#endif
{
  for (int j = 0; j < ncols; ++j) {
    const PksCols& P = Pks[j];
    int idx = choose_cp_row_cols(P, policy, rtmed_i, sn_min);
    chosen_idx[j] = idx;
    if (idx >= 0 && idx < P.n && P.has_rtmin && P.has_rtmax) {
      lo[j] = P.rtmin[idx];
      hi[j] = P.rtmax[idx];
    } else {
      chosen_idx[j] = -1;
    }
  }
}
     // ------------------------------------------------------------------------

     // build outputs per sample (sequential, R API)
     List EIC_df_list(ncols);
     List picked_cp_list(ncols);

     for (int j = 0; j < ncols; ++j) {
       S4 xChrom = get_cell_flex(input_list, i, j, nrows);
       NumericMatrix pks = as<NumericMatrix>(xChrom.slot("chromPeaks"));

       int idx = chosen_idx[j];
       if (idx < 0) {
         EIC_df_list[j] = DataFrame::create(
           Named("rt") = NumericVector(0),
           Named("intensity") = NumericVector(0)
         );
         picked_cp_list[j] = R_NilValue;
       } else {
         NumericVector rt  = xChrom.slot("rtime");
         NumericVector inx = xChrom.slot("intensity");
         EIC_df_list[j]    = subset_eic(rt, inx, lo[j], hi[j], pad);
         picked_cp_list[j] = cp_row_df(pks, idx);
       }
     }

     EIC_df_list.attr("names")    = nm_eic;
     picked_cp_list.attr("names") = nm_cp;

     out[i] = List::create(
       Named("chrom_peaks_df") = chrom_peaks_df,
       Named("EIC_df_list")    = EIC_df_list,
       Named("picked_cp_list") = picked_cp_list,
       Named("policy_used")    = policy,
       Named("rtmed")          = rtmed_i,
       Named("feature_name")   = fname
     );
   }

   return out;
 }
