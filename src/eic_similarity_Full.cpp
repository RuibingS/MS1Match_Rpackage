// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
#ifdef _OPENMP
#include <omp.h>
#endif

// ---------- small helpers ----------
static inline std::vector<double> rank_average(const std::vector<double>& x){
  const size_t n=x.size();
  std::vector<size_t> idx(n); std::iota(idx.begin(),idx.end(),0);
  std::stable_sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return x[a]<x[b];});
  std::vector<double> r(n); size_t i=0;
  while(i<n){
    size_t j=i+1; while(j<n && x[idx[j]]==x[idx[i]]) ++j;
    double rr = ((double)i+(double)j-1.0)/2.0 + 1.0;
    for(size_t k=i;k<j;++k) r[idx[k]] = rr;
    i=j;
  }
  return r;
}

static inline double corr_pearson(const std::vector<double>& a,
                                  const std::vector<double>& b){
  const size_t n=a.size();
  if(n!=b.size() || n<3) return NA_REAL;
  double ma=std::accumulate(a.begin(),a.end(),0.0)/n;
  double mb=std::accumulate(b.begin(),b.end(),0.0)/n;
  double num=0.0, da2=0.0, db2=0.0;
  for(size_t i=0;i<n;++i){
    double da=a[i]-ma, db=b[i]-mb;
    num+=da*db; da2+=da*da; db2+=db*db;
  }
  if(da2<=0.0 || db2<=0.0) return NA_REAL;
  return num/std::sqrt(da2*db2);
}

static inline double eic_corr_align(const std::vector<double>& rt1,
                                    const std::vector<double>& in1,
                                    const std::vector<double>& rt2,
                                    const std::vector<double>& in2,
                                    double rt_thr,
                                    const std::string& method){
  std::vector<double> v1, v2; v1.reserve(rt1.size()); v2.reserve(rt1.size());
  for(size_t i=0;i<rt1.size();++i){
    double y1 = in1[i];
    if(std::isnan(y1) || y1<=0) continue;
    double best = rt_thr + 1.0;
    int jbest = -1;
    for(size_t j=0;j<rt2.size();++j){
      double y2 = in2[j];
      if(std::isnan(y2) || y2<=0) continue;
      double d = std::fabs(rt1[i]-rt2[j]);
      if(d<=rt_thr && d<best){ best=d; jbest=(int)j; }
    }
    if(jbest>=0){
      v1.push_back(std::log10(y1));
      v2.push_back(std::log10(in2[(size_t)jbest]));
    }
  }
  if(v1.size()<3) return NA_REAL;
  if(method=="spearman"){ auto r1=rank_average(v1), r2=rank_average(v2); return corr_pearson(r1,r2); }
  return corr_pearson(v1,v2);
}

static inline bool is_dataframe(SEXP x){ return Rf_isObject(x) && Rf_inherits(x, "data.frame"); }

// ---------- main: DataFrame with feature_1 / feature_2 ----------
// [[Rcpp::export]]
 SEXP EIC_cor_pairs_by_sample(
     List feature_EIC_process,
     DataFrame feature_pairs,
     double rt_EIC_thr = 0.2,
     std::string method = "pearson",
     bool parallel = true,
     int n_threads = 0
 ){
   if(feature_EIC_process.size()==0) stop("feature_EIC_process is empty.");

   // --- pull pair columns safely ---
   if(!feature_pairs.containsElementNamed("feature_1") ||
      !feature_pairs.containsElementNamed("feature_2"))
      stop("feature_pairs must contain columns 'feature_1' and 'feature_2'.");

   CharacterVector f1 = as<CharacterVector>(feature_pairs["feature_1"]);
   CharacterVector f2 = as<CharacterVector>(feature_pairs["feature_2"]);
   if(f1.size()!=f2.size()) stop("'feature_1' and 'feature_2' lengths differ.");
   const int n_pairs = f1.size();

   // --- determine #samples from first usable entry ---
   int n_samp = -1;
   for(int i=0;i<feature_EIC_process.size();++i){
     List it = feature_EIC_process[i];
     if(it.containsElementNamed("EIC_df_list")){
       SEXP s_list = it["EIC_df_list"];
       if(!Rf_isNull(s_list)){
         List l = as<List>(s_list);
         n_samp = l.size();
         break;
       }
     }
   }
   if(n_samp <= 0) stop("Could not determine number of samples from EIC_df_list.");

   // --- output column names per sample ---
   std::vector<std::string> out_colnames(n_samp);
   for(int s=0;s<n_samp;++s)
     out_colnames[s] = "EIC_similarity_sample_" + std::to_string(s+1);

   struct Series { std::vector<double> rt,in; bool ok=false; };
   struct FeatData { std::vector<Series> per_sample; };
   std::unordered_map<std::string, FeatData> cache;
   cache.reserve(feature_EIC_process.size());

   CharacterVector outer_names;
   {
     SEXP nm = Rf_getAttrib(feature_EIC_process, R_NamesSymbol);
     if(Rf_isNull(nm)) outer_names = CharacterVector(0);
     else outer_names = as<CharacterVector>(nm);
   }

   for(int i=0;i<feature_EIC_process.size();++i){
     List item = feature_EIC_process[i];

     // feature key preference: inner $feature_name -> outer list name -> FT%04d
     std::string fname;
     if(item.containsElementNamed("feature_name")){
       CharacterVector fn = item["feature_name"];
       if(fn.size()>0 && fn[0]!=NA_STRING) fname = as<std::string>(fn[0]);
     }
     if(fname.empty() && outer_names.size()==feature_EIC_process.size()){
       Rcpp::String s = outer_names[i];
       if(s!=NA_STRING) fname = (std::string)s;
     }
     if(fname.empty()){ char buf[16]; std::snprintf(buf,sizeof(buf),"FT%04d", i+1); fname=std::string(buf); }

     FeatData fd; fd.per_sample.resize(n_samp);
     SEXP s_list = item.containsElementNamed("EIC_df_list") ? item["EIC_df_list"] : R_NilValue;
     if(!Rf_isNull(s_list)){
       List eic_list = as<List>(s_list);
       for(int s=0;s<n_samp;++s){
         Series ss;
         if(s<eic_list.size() && !Rf_isNull(eic_list[s]) && is_dataframe(eic_list[s])){
           DataFrame df = as<DataFrame>(eic_list[s]);
           if(df.containsElementNamed("rt") && df.containsElementNamed("intensity")){
             ss.rt = as< std::vector<double> >(df["rt"]);
             ss.in = as< std::vector<double> >(df["intensity"]);
             ss.ok = (!ss.rt.empty() && ss.rt.size()==ss.in.size());
           }
         }
         fd.per_sample[s] = std::move(ss);
       }
     }
     cache.emplace(fname, std::move(fd));   // last wins if duplicate key
   }

   std::vector<std::string> pair_f1(n_pairs), pair_f2(n_pairs);
   for(int k=0;k<n_pairs;++k){
     pair_f1[k] = (f1[k]==NA_STRING) ? std::string() : (std::string)f1[k];
     pair_f2[k] = (f2[k]==NA_STRING) ? std::string() : (std::string)f2[k];
   }

   // --- outputs ---
   std::vector<std::string> out_f1 = pair_f1, out_f2 = pair_f2;
   std::vector< std::vector<double> > out_sim(n_samp, std::vector<double>(n_pairs, NA_REAL));

#ifdef _OPENMP
   if(parallel && n_threads>0) omp_set_num_threads(n_threads);
#endif

#ifdef _OPENMP
   if(parallel){
#pragma omp parallel for schedule(dynamic)
     for(int k=0;k<n_pairs;++k){
       auto ita = cache.find(pair_f1[k]);
       auto itb = cache.find(pair_f2[k]);
       if(ita!=cache.end() && itb!=cache.end()){
         for(int s=0;s<n_samp;++s){
           const Series& A = ita->second.per_sample[s];
           const Series& B = itb->second.per_sample[s];
           if(A.ok && B.ok)
             out_sim[s][k] = eic_corr_align(A.rt,A.in,B.rt,B.in,rt_EIC_thr,method);
         }
       }
     }
   } else
#endif
{
  for(int k=0;k<n_pairs;++k){
    auto ita = cache.find(pair_f1[k]);
    auto itb = cache.find(pair_f2[k]);
    if(ita!=cache.end() && itb!=cache.end()){
      for(int s=0;s<n_samp;++s){
        const Series& A = ita->second.per_sample[s];
        const Series& B = itb->second.per_sample[s];
        if(A.ok && B.ok)
          out_sim[s][k] = eic_corr_align(A.rt,A.in,B.rt,B.in,rt_EIC_thr,method);
      }
    }
  }
}

   List res(2 + n_samp);
   res[0] = out_f1;
   res[1] = out_f2;
   for(int s=0;s<n_samp;++s)
     res[2+s] = NumericVector(out_sim[s].begin(), out_sim[s].end());

   CharacterVector cn(2 + n_samp);
   cn[0]="feature_1"; cn[1]="feature_2";
   for(int s=0;s<n_samp;++s) cn[2+s] = out_colnames[s];
   res.attr("names") = cn;
   res.attr("class") = "data.frame";
   res.attr("row.names") = IntegerVector::create(NA_INTEGER, -n_pairs);
   return res;
 }
