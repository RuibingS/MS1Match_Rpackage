#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame summarize_peak_mz_cpp(
    List lib,
    double ppm_tol = 5.0,          // ppm tolerance for equality
    double isotope_shift = 1.0034, // Da difference per isotope
    int max_isotope = 3
) {
  const double sodium = 22.989218;
  const double ppm_factor = 1e-6;

  std::vector<std::string> names;
  std::vector<double> exactMasses;
  std::vector<int> n_leq, n_gt, n_total;

  for (int i = 0; i < lib.size(); ++i) {
    Rcpp::List compound = lib[i];
    if (!compound.containsElementNamed("Peaks") ||
        !compound.containsElementNamed("ExactMass"))
        continue;

    double exactMass = as<double>(compound["ExactMass"]);
    DataFrame peaks = as<DataFrame>(compound["Peaks"]);
    NumericVector mz = peaks["mz"];
    int total = mz.size();

    // --- Base reference: [M+Na]+
    double base_mz = exactMass + sodium;
    double tol_abs = base_mz * ppm_tol * ppm_factor;

    // --- Include isotopes (+1.0034 × 1..max_isotope)
    std::vector<double> isotope_mz(max_isotope);
    for (int k = 0; k < max_isotope; ++k)
      isotope_mz[k] = base_mz + (k + 1) * isotope_shift;

    int c_leq = 0, c_gt = 0;

    for (int j = 0; j < total; ++j) {
      double mz_val = mz[j];
      double diff = mz_val - base_mz;
      double ppm_diff = std::fabs(diff) / base_mz * 1e6;

      bool is_isotope = false;
      for (int k = 0; k < max_isotope; ++k) {
        if (std::fabs(mz_val - isotope_mz[k]) <= tol_abs) {
          is_isotope = true;
          break;
        }
      }

      if (mz_val <= base_mz + tol_abs || is_isotope) {
        c_leq++;
      } else {
        c_gt++;
      }
    }

    std::string name = "";
    if (compound.containsElementNamed("Name"))
      name = as<std::string>(compound["Name"]);

    names.push_back(name);
    exactMasses.push_back(exactMass);
    n_leq.push_back(c_leq);
    n_gt.push_back(c_gt);
    n_total.push_back(total);
  }

  return DataFrame::create(
    Named("Compound") = names,
    Named("ExactMass") = exactMasses,
    Named("Adduct") = "[M+Na]+",
    Named("n_smaller_M_Na") = n_leq,
    Named("n_bigger_M_Na") = n_gt,
    Named("n_total") = n_total
  );
}
