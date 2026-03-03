#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List compare_library_peaks(List in_house_list,
                                     List mona_lib_list,
                                     DataFrame overlaps_df,
                                     double mz_tol = 0.005) {

  CharacterVector names_inhouse = overlaps_df["Name_in_house"];
  CharacterVector names_mona    = overlaps_df["Name_Massbank"];
  int n = names_inhouse.size();

  std::vector<std::string> name_inhouse_out;
  std::vector<std::string> name_mona_out;
  std::vector<int> n_inhouse_only;
  std::vector<int> n_mona_only;
  std::vector<int> n_overlap;
  List mz_inhouse_only_list;
  List mz_mona_only_list;
  List mz_overlap_list;

  auto get_mz = [](List entry) -> NumericVector {
    if (!entry.containsElementNamed("Peaks")) return NumericVector();
    DataFrame df = as<DataFrame>(entry["Peaks"]);
    if (!df.containsElementNamed("mz")) return NumericVector();
    return df["mz"];
  };


  for (int k = 0; k < n; ++k) {
    std::string name_inhouse = as<std::string>(names_inhouse[k]);
    std::string name_mona    = as<std::string>(names_mona[k]);

    // Find in-house entry
    List in_entry;
    bool found_in = false;
    for (int i = 0; i < in_house_list.size(); ++i) {
      List tmp = in_house_list[i];
      if (tmp.containsElementNamed("Name")) {
        std::string n_tmp = as<std::string>(tmp["Name"]);
        if (n_tmp == name_inhouse) {
          in_entry = tmp;
          found_in = true;
          break;
        }
      }
    }
    if (!found_in) continue;

    // Find MoNA entry
    List mona_entry;
    bool found_mona = false;
    for (int j = 0; j < mona_lib_list.size(); ++j) {
      List tmp = mona_lib_list[j];
      if (tmp.containsElementNamed("Name")) {
        std::string n_tmp = as<std::string>(tmp["Name"]);
        if (n_tmp == name_mona) {
          mona_entry = tmp;
          found_mona = true;
          break;
        }
      }
    }
    if (!found_mona) continue;

    NumericVector mz_in = get_mz(in_entry);
    NumericVector mz_mona = get_mz(mona_entry);

    NumericVector overlap_vals, in_unique_vals, mona_unique_vals;
    if (mz_in.size() > 0 && mz_mona.size() > 0) {
      LogicalVector matched_mona(mz_mona.size(), false);
      for (int i = 0; i < mz_in.size(); ++i) {
        double val = mz_in[i];
        bool matched = false;
        for (int j = 0; j < mz_mona.size(); ++j) {
          if (matched_mona[j]) continue;
          if (std::fabs(val - mz_mona[j]) <= mz_tol) {
            matched = true;
            matched_mona[j] = true;
            overlap_vals.push_back(val);
            break;
          }
        }
        if (!matched) in_unique_vals.push_back(val);
      }
      for (int j = 0; j < mz_mona.size(); ++j)
        if (!matched_mona[j]) mona_unique_vals.push_back(mz_mona[j]);
    } else {
      in_unique_vals = mz_in;
      mona_unique_vals = mz_mona;
    }

    name_inhouse_out.push_back(name_inhouse);
    name_mona_out.push_back(name_mona);
    n_inhouse_only.push_back(in_unique_vals.size());
    n_mona_only.push_back(mona_unique_vals.size());
    n_overlap.push_back(overlap_vals.size());
    mz_inhouse_only_list.push_back(in_unique_vals);
    mz_mona_only_list.push_back(mona_unique_vals);
    mz_overlap_list.push_back(overlap_vals);
  }

  return List::create(
    Named("Name_in_house")  = name_inhouse_out,
    Named("Name_Massbank")  = name_mona_out,
    Named("Unique_inhouse") = n_inhouse_only,
    Named("Unique_MoNA")    = n_mona_only,
    Named("Overlap_peaks")  = n_overlap,
    Named("mz_inhouse_only")= mz_inhouse_only_list,
    Named("mz_mona_only")   = mz_mona_only_list,
    Named("mz_overlap")     = mz_overlap_list
  );
}
