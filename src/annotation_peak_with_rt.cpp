#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iomanip>
#include "utils_check.h"
#include "utils_sort.h"
#include <regex>
#include <sstream>
#include "utils_annotation.h"


using namespace Rcpp;


//prototypes

Rcpp::List PImatch_fun(Rcpp::DataFrame FT,
                       Rcpp::DataFrame Comp_data,
                       Rcpp::CharacterVector PIon,
                       double diff_mz_thr,
                       double diff_rt_thr,
                       int ppm);

Rcpp::List findNL_fun(
    Rcpp::DataFrame NL_data,
    Rcpp::DataFrame adj_long,
    double diff_mz_thr,
    double diff_rt_thr,
    int ppm);



CharacterVector Check_Iso_Charge(CharacterVector name,
                                 NumericVector mz,
                                 NumericVector rt,
                                 NumericVector intensity,
                                 double iso_threshold ,
                                 int max_charge ,
                                 int max_isotope);


double compute_adduct_mass(double M, std::string adduct_formula);


// [[Rcpp::export]]
List annotate_with_rt_cpp(List compounds_list, Rcpp::DataFrame FT,
                           Rcpp::DataFrame adducts_df,
                           Rcpp::DataFrame neutral_loss_table,
                           std::string polarity,
                           NumericVector intensity,
                           double mz_threshold = 0.005, double PI_rt_threshold = 5.0,
                           double rt_threshold = 3.0,
                           double iso_threshold = 1.003355, int max_charge = 3, int max_isotope = 4,
                           bool run_PI_match = true, bool run_adduct_match = true, bool run_NL_match = true) {

  NumericVector FT_mz = get_column_with_fallback<NumericVector>(FT, "mzmed", "mz", "m/z");
  CharacterVector feature_name = get_required_column<CharacterVector>(FT, "feature_name");
  NumericVector FT_rt = get_column_with_fallback<NumericVector>(FT, "rtmed", "rt", "retention time");





  bool is_sorted = true;
  for (int i = 1; i < FT_mz.size(); ++i) {
    if (FT_mz[i] < FT_mz[i - 1]) {
      is_sorted = false;
      break;
    }
  }
  if (!is_sorted) {
    Rcpp::Rcout << "Sorting FT by m/z...\n";
    std::vector<int> idx = compute_sorted_index(FT_mz);
    FT_mz = reorder_vector(FT_mz, idx);
    FT_rt = reorder_vector(FT_rt, idx);
    feature_name = reorder_vector(feature_name, idx);
  }

  int n_features = FT.nrows();
  int n_compounds = compounds_list.size();
  Rcpp::Rcout << "Features: " << n_features << ", Compounds: " << n_compounds << "\n";

  std::vector<double> compound_mass(n_compounds);
  std::vector<double> compound_rt(n_compounds);
  std::vector<std::string> compound_name(n_compounds);
  for (int i = 0; i < n_compounds; ++i) {
    List compound = compounds_list[i];
    if (!compound.containsElementNamed("ExactMass") || !compound.containsElementNamed("RetTime")) continue;
    compound_mass[i] = compound["ExactMass"];
    compound_rt[i] = std::atof(as<std::string>(compound["RetTime"]).c_str());
    compound_name[i] = as<std::string>(compound["Name"]);
  }

  CharacterVector adduct_ions = adducts_df["Ion.name"];
  CharacterVector adduct_polarities = adducts_df["Ion.polarity"];
  CharacterVector adduct_formulas = adducts_df["Ion.mass"];

  NumericVector nl_mz_difference = neutral_loss_table["Accurate.Mass"];
  CharacterVector nl_name = neutral_loss_table["Formula"];
  CharacterVector pos_polarity = neutral_loss_table["Pos"];
  CharacterVector neg_polarity = neutral_loss_table["Neg"];

  std::vector<std::string> PI_annotations(n_features);
  std::vector<std::string> adduct_annotations(n_features);
  std::vector<std::string> NL_annotations(n_features);

  if (run_PI_match) {
    Rcpp::Rcout << "Starting PI match...\n";
    for (int f = 0; f < n_features; ++f) {
      double mz = FT_mz[f];
      double rt = FT_rt[f];

      for (int i = 0; i < n_compounds; ++i) {
        double M = compound_mass[i];
        double crt = compound_rt[i];
        if (std::abs(rt - crt) > PI_rt_threshold) continue;

        if (polarity == "pos") {
          if (std::abs(mz - (M + 1.007276)) <= mz_threshold) {
            PI_annotations[f] = compound_name[i] + " [M+H]+";
            break;
          }
        } else if (polarity == "neg") {
          if (std::abs(mz - (M - 1.007276)) <= mz_threshold) {
            PI_annotations[f] = compound_name[i] + " [M-H]-";
            break;
          }
        }
      }
    }
  }
  if (run_adduct_match) {
    Rcpp::Rcout << "Starting adduct match...\n";

    //  Filtered adduct table by ion polarity
    std::vector<std::string> filtered_formulas;
    std::vector<std::string> filtered_ions;
    for (int i = 0; i < adduct_formulas.size(); ++i) {
      if (adduct_polarities[i] == polarity) {
        filtered_formulas.push_back(Rcpp::as<std::string>(adduct_formulas[i]));
        filtered_ions.push_back(Rcpp::as<std::string>(adduct_ions[i]));
      }
    }

    // For each feature
    for (int f = 0; f < n_features; ++f) {
      double mz_feature = FT_mz[f];
      std::vector<std::string> matches;

      // For each compound
      for (int c = 0; c < n_compounds; ++c) {
        double M = compound_mass[c];
        std::string cname = compound_name[c];

        // Calculate theoretical adduct m/z values for each M
        std::vector<double> mz_adducts = parse_adduct_vector(M, filtered_formulas);

        for (size_t k = 0; k < mz_adducts.size(); ++k) {
          double mz_adduct = mz_adducts[k];
          double diff = std::abs(mz_feature - mz_adduct);
          double mz_tol = mz_threshold;

          if (diff <= mz_tol) {
            std::string annotation = filtered_ions[k] + " " + std::to_string(M); // cname + " " + filtered_ions[k] + " " + std::to_string(M);
            matches.push_back(annotation);
          }
        }
      }

      // Concatenate  all adduct matches into a string
      adduct_annotations[f] = join_annotations(matches);
    }

    Rcpp::Rcout << "Adduct annotation completed.\n";
  }



  if (run_NL_match) {
    Rcpp::Rcout << "Starting neutral loss match...\n";
    for (int f1 = 0; f1 < n_features; ++f1) {
      double mz1 = FT_mz[f1];
      double rt1 = FT_rt[f1];
      std::vector<std::string> nl_matches;

      for (int f2 = 0; f2 < n_features; ++f2) {
        if (f1 == f2) continue;
        double mz2 = FT_mz[f2];
        double rt2 = FT_rt[f2];
        if (std::abs(rt1 - rt2) > rt_threshold || mz2 <= mz1) continue;

        double diff = mz2 - mz1;
        for (int m = 0; m < nl_mz_difference.size(); ++m) {
          bool match = ((polarity == "pos" && as<std::string>(pos_polarity[m]) == "+") ||
                        (polarity == "neg" && as<std::string>(neg_polarity[m]) == "+"));
          if (match && std::abs(diff - nl_mz_difference[m]) <= mz_threshold) {
            nl_matches.push_back("-" + as<std::string>(nl_name[m]) + " (" + as<std::string>(feature_name[f2]) + ")");
          }
        }
      }
      NL_annotations[f1] = join_annotations(nl_matches);
    }
  }
  Rcpp::Rcout << "Starting isotope anotation...\n";
  CharacterVector isotope_labels = Check_Iso_Charge(feature_name, FT_mz, FT_rt, intensity, iso_threshold, max_charge, max_isotope);
  Rcpp::Rcout << "Isotope anotation completed\n";

  List output_list(n_features);
  for (int i = 0; i < n_features; ++i) {
    output_list[i] = List::create(
      Named("feature_name") = CharacterVector::create(feature_name[i]),
      Named("mz") = FT_mz[i],
      Named("rt") = FT_rt[i],
      Named("PI_match") = CharacterVector::create(PI_annotations[i]),
      Named("adduct_match") = CharacterVector::create(adduct_annotations[i]),
      Named("NL_match") = NL_annotations[i],
      Named("iso_annotation") = CharacterVector::create(isotope_labels[i])
    );


  }

  Rcpp::Rcout << "Annotation completed.\n";
  return output_list;
}









// Isotope annotation

// [[Rcpp::export]]
Rcpp::CharacterVector Check_Iso_Charge(Rcpp::CharacterVector name,
                                       Rcpp::NumericVector mz,
                                       Rcpp::NumericVector rt,
                                       Rcpp::NumericVector intensity,
                                       double iso_threshold = 1.003355,
                                       int max_charge = 3,
                                       int max_isotope = 4) {
  // Just warn if mz is not sorted
  bool is_sorted = true;
  for (int i = 1; i < mz.size(); ++i) {
    if (mz[i] < mz[i - 1]) {
      is_sorted = false;
      break;
    }
  }

  if (!is_sorted) {
    Rcpp::Rcout << "Warning: input m/z is not sorted in ascending order. Please make sure the feature table is sorted by m/z.\n";
  }

  Rcpp::CharacterVector labels(name.size());
  int groupIndex = 1;

  for (int i = 0; i < name.size(); i++) {
    if (labels[i] != "") continue;

    std::vector<int> maxIsotopeForCharge(max_charge + 1, 0);

    for (int j = i + 1; j < name.size(); j++) {
      if (labels[j] != "") continue;

      double mz_diff = mz[j] - mz[i];
      double rt_diff = std::abs(rt[j] - rt[i]);

      for (int c = 1; c <= max_charge; c++) {
        for (int m = 1; m <= max_isotope; m++) {
          double expected_mz_diff = m * iso_threshold / c;
          double max_int_ratio = (mz[i] * c / 12) * 0.011;

          if (m == 1) {
            if (std::abs(mz_diff - expected_mz_diff) <= 0.002 &&
                rt_diff < 3 &&
                intensity[i] * max_int_ratio >= intensity[j]) {

              labels[i] = Rcpp::String("[" + std::to_string(groupIndex) + "][M]" +
                (c > 1 ? std::to_string(c) : "") + "+");

              labels[j] = Rcpp::String("[" + std::to_string(groupIndex) + "][M+" +
                std::to_string(m) + "]" +
                (c > 1 ? std::to_string(c) : "") + "+");

              maxIsotopeForCharge[c] = m;
            }
          } else {
            if (std::abs(mz_diff - expected_mz_diff) <= 0.01 &&
                rt_diff < 3 &&
                intensity[j] <= intensity[i]) {
              if (maxIsotopeForCharge[c] == m - 1) {
                labels[j] = Rcpp::String("[" + std::to_string(groupIndex) + "][M+" +
                  std::to_string(m) + "]" +
                  (c > 1 ? std::to_string(c) : "") + "+");

                maxIsotopeForCharge[c] = m;
              }
            }
          }
        }
      }
    }

    if (!labels[i].empty()) {
      groupIndex++;
    }
  }

  return labels;
}




// [[Rcpp::export]]
Rcpp::List PImatch_fun(Rcpp::DataFrame FT,
                       Rcpp::DataFrame Comp_data,
                       Rcpp::CharacterVector PIon,
                       double diff_mz_thr,
                       double diff_rt_thr,
                       int ppm) {
  int nrow_FT = FT.nrow();
  int nrow_Comp = Comp_data.nrow();
  int length_PIon = PIon.length();

  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mzmed"];
  Rcpp::NumericVector FT_rt = FT["rtmed"];

  Rcpp::CharacterVector Comp_name = Comp_data["Name"];
  Rcpp::NumericVector Comp_rt = Comp_data["RT"];
  Rcpp::NumericVector Comp_M = Comp_data["exactMass"];

  List output(nrow_FT);
  for(int i = 0; i < nrow_FT ; i++)
  {
    Rcpp::List df(6);
    df.names() = Rcpp::CharacterVector::create("Feature_name", "Feature_mz","Feature_rt","Comp_name","PI_name","M");

    CharacterVector feaNames;
    CharacterVector CompNames;
    CharacterVector PIname;
    Rcpp::NumericVector  fea_rt;
    Rcpp::NumericVector  fea_mz;
    Rcpp::NumericVector  M_vec;

    for(int p=0; p<length_PIon ; p++){
      String PIon_name = PIon[p];
      Rcpp::NumericVector comp_mz = Comp_data[PIon_name];

      for(int j = 0; j < nrow_Comp; j++){

        float diff_mz= std::fabs(FT_mz[i]-comp_mz[j]);

        if (((Comp_rt[j] != 0 || !NumericVector::is_na(Comp_rt[j])) &&
            (diff_mz <= comp_mz[j] * ppm * 1e-6) &&
            (std::abs(FT_rt[i] - Comp_rt[j]) <= diff_rt_thr)) ||
            ((diff_mz <= diff_mz_thr) &&
            (std::abs(FT_rt[i] - Comp_rt[j]) <= diff_rt_thr)))
        {
          feaNames.push_back(feature_name[i]);
          fea_rt.push_back(FT_rt[i]);
          fea_mz.push_back(FT_mz[i]);
          CompNames.push_back(Comp_name[j]);
          PIname.push_back(PIon[p]);
          M_vec.push_back(Comp_M[j]);


        }

        else{
          if(diff_mz<= comp_mz[j]*ppm* 1e-6||diff_mz<= diff_mz_thr )
          {
            feaNames.push_back(feature_name[i]);
            fea_rt.push_back(FT_rt[i]);
            fea_mz.push_back(FT_mz[i]);
            CompNames.push_back(Comp_name[j]);
            PIname.push_back(PIon[p]);
            M_vec.push_back(Comp_M[j]);

          }
        }

      }

    }

    df["Feature_name"] = feaNames;
    df["Feature_mz"] = fea_mz;
    df["Feature_rt"] = fea_rt;
    df["Comp_name"] = CompNames;
    df["PI_name"] = PIname;
    df["M"] = M_vec;
    output[i] = df;
  }
  return output;

}



// NL match
// [[Rcpp::export]]
Rcpp::List findNL_fun(
    Rcpp::DataFrame NL_data,
    Rcpp::DataFrame adj_long,
    double diff_mz_thr,
    double diff_rt_thr,
    int ppm) {

  int nrow_adj_long = adj_long.nrow();
  int nrow_NL = NL_data.nrow();

  Rcpp::NumericVector mz_diff = adj_long["mz_diff"];
  Rcpp::NumericVector rt_diff = adj_long["rt_diff"];
  Rcpp::NumericVector NL_mz = NL_data["Accurate.Mass"];
  Rcpp::CharacterVector NL_name = NL_data["Neutral.Loss"];
  Rcpp::CharacterVector fea_name1 = adj_long["Var1"];
  Rcpp::CharacterVector fea_name2 = adj_long["Var2"];

  Rcpp::NumericVector mz1 = adj_long["mz_x"];
  Rcpp::NumericVector mz2 = adj_long["mz_y"];

  List out(nrow_adj_long);

  for (int i = 0; i < nrow_adj_long - 1; i++) {
    Rcpp::List df(3);
    df.names() = Rcpp::CharacterVector::create("Feature_name1", "Feature_name2", "NL_Names");
    CharacterVector feaNames1;
    CharacterVector feaNames2;
    CharacterVector NLNames;

    for (int j = 0; j < nrow_NL - 1; j++) {
      float diff_mz1 = std::fabs(NL_mz[j] - mz_diff[i]);

      if ((diff_mz1 <= diff_mz_thr || diff_mz1 <= NL_mz[j] * ppm * 1e-6) && rt_diff[i] <= diff_rt_thr) {


        if (mz1[i] < mz2[i]) {
          feaNames1.push_back(fea_name2[i]); // swap names
          feaNames2.push_back(fea_name1[i]);
        } else {
          feaNames1.push_back(fea_name1[i]);
          feaNames2.push_back(fea_name2[i]);
        }

        NLNames.push_back(NL_name[j]);
      }
    }

    df["Feature_name1"] = feaNames1;
    df["Feature_name2"] = feaNames2;
    df["NL_Names"] = NLNames;
    out[i] = df;
  }

  return out;
}


