#include <Rcpp.h>
#include <unordered_set>
#include <cmath>
#include <sstream>
using namespace Rcpp;

// helper: bin m/z to avoid duplicates due to tolerance
inline double round_mz(double val, double tol) {
  return std::round(val / tol) * tol;
}

// [[Rcpp::export]]
List compare_library_peaks_3lib(List in_house_list,
                                List mona_lib_list,
                                List csu_lib_list,
                                DataFrame overlaps_df,
                                double mz_tol = 0.005) {

  CharacterVector names_inhouse = overlaps_df["Name_CBIO"];
  CharacterVector names_mona    = overlaps_df["Name_MassBank"];
  CharacterVector names_csu     = overlaps_df["Name_CSU"];
  int n = names_inhouse.size();

  // output containers
  std::vector<std::string> name_inhouse_out, name_mona_out, name_csu_out;
  std::vector<int> n_unique_in, n_unique_mo, n_unique_csu;
  std::vector<int> n_overlap_in_mo, n_overlap_in_csu, n_overlap_mo_csu, n_overlap_3way;

  List mz_in_only_list, mz_mo_only_list, mz_csu_only_list;
  List mz_in_mo_list, mz_in_csu_list, mz_mo_csu_list, mz_3way_list;

  // helper: extract m/z vector
  auto get_mz = [](List entry) -> NumericVector {
    if (!entry.containsElementNamed("Peaks")) return NumericVector();
    DataFrame df = as<DataFrame>(entry["Peaks"]);
    if (!df.containsElementNamed("mz")) return NumericVector();
    return df["mz"];
  };

  for (int k = 0; k < n; ++k) {
    std::string name_in  = as<std::string>(names_inhouse[k]);
    std::string name_mo  = as<std::string>(names_mona[k]);
    std::string name_csu = as<std::string>(names_csu[k]);

    // find entries
    List in_entry, mo_entry, csu_entry;
    bool found_in = false, found_mo = false, found_csu = false;

    for (int i = 0; i < in_house_list.size(); ++i) {
      List tmp = in_house_list[i];
      if (tmp.containsElementNamed("Name") &&
          as<std::string>(tmp["Name"]) == name_in) {
        in_entry = tmp; found_in = true; break;
      }
    }

    for (int j = 0; j < mona_lib_list.size(); ++j) {
      List tmp = mona_lib_list[j];
      if (tmp.containsElementNamed("Name") &&
          as<std::string>(tmp["Name"]) == name_mo) {
        mo_entry = tmp; found_mo = true; break;
      }
    }

    for (int h = 0; h < csu_lib_list.size(); ++h) {
      List tmp = csu_lib_list[h];
      if (tmp.containsElementNamed("Name") &&
          as<std::string>(tmp["Name"]) == name_csu) {
        csu_entry = tmp; found_csu = true; break;
      }
    }

    if (!found_in || !found_mo || !found_csu) {
      Rcout << "Skipping compound (not found in all libs): "
            << name_in << std::endl;
      continue;
    }

    NumericVector mz_in  = get_mz(in_entry);
    NumericVector mz_mo  = get_mz(mo_entry);
    NumericVector mz_csu = get_mz(csu_entry);

    if (mz_csu.size() > 500) {
      Rcout << "Warning: " << name_csu
            << " in CSU_PMF has " << mz_csu.size()
            << " peaks!" << std::endl;
    }

    // containers
    std::unordered_set<double> seen_in, seen_mo, seen_csu;
    NumericVector in_only, mo_only, csu_only;
    NumericVector in_mo, in_csu, mo_csu, all3;

    // ---- check Inhouse peaks ----
    for (double val_in : mz_in) {
      bool match_m = false, match_c = false;

      for (double val_mo : mz_mo)
        if (std::fabs(val_in - val_mo) <= mz_tol) { match_m = true; break; }

        for (double val_csu : mz_csu)
          if (std::fabs(val_in - val_csu) <= mz_tol) { match_c = true; break; }

          double key = round_mz(val_in, mz_tol);
          if (match_m && match_c)
            all3.push_back(val_in);
          else if (match_m)
            in_mo.push_back(val_in);
          else if (match_c)
            in_csu.push_back(val_in);
          else
            in_only.push_back(val_in);
          seen_in.insert(key);
    }

    // ---- check MoNA peaks ----
    for (double val_mo : mz_mo) {
      double key = round_mz(val_mo, mz_tol);
      if (seen_mo.count(key)) continue;
      seen_mo.insert(key);

      bool in_in = false, in_c = false;

      for (double val_in : mz_in)
        if (std::fabs(val_mo - val_in) <= mz_tol) { in_in = true; break; }

        for (double val_csu : mz_csu)
          if (std::fabs(val_mo - val_csu) <= mz_tol) { in_c = true; break; }

          if (in_in && in_c)
            continue; // already counted as 3-way
          else if (in_in)
            in_mo.push_back(val_mo);
          else if (in_c)
            mo_csu.push_back(val_mo);
          else
            mo_only.push_back(val_mo);
    }

    // ---- check CSU peaks ----
    for (double val_csu : mz_csu) {
      double key = round_mz(val_csu, mz_tol);
      if (seen_csu.count(key)) continue;
      seen_csu.insert(key);

      bool in_in = false, in_m = false;

      for (double val_in : mz_in)
        if (std::fabs(val_csu - val_in) <= mz_tol) { in_in = true; break; }

        for (double val_mo : mz_mo)
          if (std::fabs(val_csu - val_mo) <= mz_tol) { in_m = true; break; }

          if (in_in && in_m)
            continue; // 3-way
          else if (in_in)
            in_csu.push_back(val_csu);
          else if (in_m)
            mo_csu.push_back(val_csu);
          else
            csu_only.push_back(val_csu);
    }

    // ---- Debug summary ----
    Rcout << "Compound: " << name_in
          << " | Peaks Inhouse: " << mz_in.size()
          << ", MoNA: " << mz_mo.size()
          << ", CSU_PMF: " << mz_csu.size()
          << " | Unique_in: " << in_only.size()
          << ", Unique_mo: " << mo_only.size()
          << ", Unique_csu: " << csu_only.size()
          << ", in∩mo: " << in_mo.size()
          << ", in∩csu: " << in_csu.size()
          << ", mo∩csu: " << mo_csu.size()
          << ", 3-way: " << all3.size() << std::endl;

    // ---- store results ----
    name_inhouse_out.push_back(name_in);
    name_mona_out.push_back(name_mo);
    name_csu_out.push_back(name_csu);

    n_unique_in.push_back(in_only.size());
    n_unique_mo.push_back(mo_only.size());
    n_unique_csu.push_back(csu_only.size());
    n_overlap_in_mo.push_back(in_mo.size());
    n_overlap_in_csu.push_back(in_csu.size());
    n_overlap_mo_csu.push_back(mo_csu.size());
    n_overlap_3way.push_back(all3.size());

    mz_in_only_list.push_back(in_only);
    mz_mo_only_list.push_back(mo_only);
    mz_csu_only_list.push_back(csu_only);
    mz_in_mo_list.push_back(in_mo);
    mz_in_csu_list.push_back(in_csu);
    mz_mo_csu_list.push_back(mo_csu);
    mz_3way_list.push_back(all3);
  }

  return List::create(
    Named("Name_Inhouse") = name_inhouse_out,
    Named("Name_MoNA") = name_mona_out,
    Named("Name_CSU_PMF") = name_csu_out,
    Named("Unique_Inhouse") = n_unique_in,
    Named("Unique_MoNA") = n_unique_mo,
    Named("Unique_CSU_PMF") = n_unique_csu,
    Named("Overlap_in_mo") = n_overlap_in_mo,
    Named("Overlap_in_csu") = n_overlap_in_csu,
    Named("Overlap_mo_csu") = n_overlap_mo_csu,
    Named("Overlap_3way") = n_overlap_3way,
    Named("mz_inhouse_only") = mz_in_only_list,
    Named("mz_mona_only") = mz_mo_only_list,
    Named("mz_csu_only") = mz_csu_only_list,
    Named("mz_overlap_in_mo") = mz_in_mo_list,
    Named("mz_overlap_in_csu") = mz_in_csu_list,
    Named("mz_overlap_mo_csu") = mz_mo_csu_list,
    Named("mz_overlap_3way") = mz_3way_list
  );
}
