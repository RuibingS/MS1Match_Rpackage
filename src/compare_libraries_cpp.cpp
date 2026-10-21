#include <Rcpp.h>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <sstream>
#include <algorithm>
using namespace Rcpp;

// =============================================================================
// compare_libraries_cpp.cpp
//
// Fast Rcpp backends for spectral library comparison.
//
// Exported functions:
//   classify_match_cpp()      — vectorised match classification (used by
//                               compare_libraries_3 to replace apply() loop)
//   compare_libraries_cpp()   — full pairwise comparison with hash-map on
//                               Formula (replaces nested R loop)
//
// All string normalisation (lowercase, trim, UTF-8) is done on the R side
// before calling these functions. C++ only does exact std::string comparison.
// =============================================================================


// -----------------------------------------------------------------------------
// Internal helpers — not exported to R
// -----------------------------------------------------------------------------

// Split a ";"-separated string into a set of trimmed non-empty tokens
static std::set<std::string> split_to_set(const std::string& s) {
  std::set<std::string> out;
  if (s.empty()) return out;
  std::stringstream ss(s);
  std::string token;
  while (std::getline(ss, token, ';')) {
    size_t start = token.find_first_not_of(" \t");
    size_t end   = token.find_last_not_of(" \t");
    if (start != std::string::npos)
      out.insert(token.substr(start, end - start + 1));
  }
  return out;
}

// Return true if two sets share at least one element
static bool sets_intersect(const std::set<std::string>& a,
                           const std::set<std::string>& b) {
  for (const auto& x : a)
    if (b.count(x)) return true;
    return false;
}

// Validate InChIKey format: 27 chars, XXXXXXXXXXXXXX-XXXXXXXXXX-X (uppercase)
static bool valid_inchikey(const std::string& ik) {
  if (ik.size() != 27)      return false;
  if (ik[14] != '-')        return false;
  if (ik[25] != '-')        return false;
  for (int i = 0;  i < 14; ++i) if (!std::isupper((unsigned char)ik[i])) return false;
  for (int i = 15; i < 25; ++i) if (!std::isupper((unsigned char)ik[i])) return false;
  if (!std::isupper((unsigned char)ik[26])) return false;
  return true;
}

// Extract first 14-char connectivity block from a validated InChIKey
static std::string ik_block1(const std::string& ik) {
  return ik.substr(0, 14);
}

// Classify a single compound pair given pre-normalised strings.
// Inputs:  ik1/ik2    — InChIKey (uppercase, "" if absent)
//          name1/name2 — Name (lowercased + trimmed)
//          syn1/syn2   — ";"-separated synonyms (lowercased + trimmed, "" if none)
// Returns: match type string
static std::string classify_pair(
    const std::string& ik1,   const std::string& ik2,
    const std::string& name1, const std::string& name2,
    const std::string& syn1,  const std::string& syn2
) {
  bool ik1v = valid_inchikey(ik1);
  bool ik2v = valid_inchikey(ik2);

  // Step 2a: full InChIKey match
  if (ik1v && ik2v && ik1 == ik2)
    return "Formula+InChIKey";

  // Step 2b: InChIKey connectivity block (first 14 chars)
  if (ik1v && ik2v && ik_block1(ik1) == ik_block1(ik2))
    return "Formula+InChIKey_connectivity";

  // Step 2c: exact Name match (already lowercased on R side)
  if (!name1.empty() && !name2.empty() && name1 == name2)
    return "Formula+Name";

  // Step 2d: synonym match — build full name+synonym sets and intersect
  std::set<std::string> set1 = split_to_set(syn1);
  std::set<std::string> set2 = split_to_set(syn2);
  if (!name1.empty()) set1.insert(name1);
  if (!name2.empty()) set2.insert(name2);

  if ((set1.size() > 1 || set2.size() > 1) && sets_intersect(set1, set2))
    return "Formula+Synonym";

  return "Formula_only";
}


// -----------------------------------------------------------------------------
// classify_match_cpp
//
// Vectorised match classification over pre-merged rows.
// Used by compare_libraries_3() to replace the slow R apply() loop.
//
// All input vectors must be the same length (one element per merged row).
// Strings must already be normalised on the R side:
//   ik1, ik2    — uppercase InChIKey or "" if absent
//   name1/2     — lowercased + trimmed name or ""
//   syn1/2      — ";"-separated lowercased synonyms or ""
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
CharacterVector classify_match_cpp(
    CharacterVector ik1,   CharacterVector ik2,
    CharacterVector name1, CharacterVector name2,
    CharacterVector syn1,  CharacterVector syn2
) {
  int n = ik1.size();
  CharacterVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = classify_pair(
      Rcpp::as<std::string>(ik1[i]),
      Rcpp::as<std::string>(ik2[i]),
      Rcpp::as<std::string>(name1[i]),
      Rcpp::as<std::string>(name2[i]),
      Rcpp::as<std::string>(syn1[i]),
      Rcpp::as<std::string>(syn2[i])
    );
  }
  return out;
}


// -----------------------------------------------------------------------------
// compare_libraries_cpp
//
// Full pairwise library comparison using std::unordered_map keyed on Formula.
// Avoids O(n*m) nested loop — only compound pairs sharing a formula are
// compared. For 600 x 700 compounds this reduces ~420,000 comparisons down
// to only the formula-matched subset.
//
// Normalised inputs (for matching — lowercased/trimmed on R side):
//   formula1, formula2  — molecular formula strings
//   ik1, ik2            — InChIKey (uppercase, "" if missing)
//   name1, name2        — Name (lowercased + trimmed, "" if missing)
//   syn1, syn2          — ";"-separated synonyms (lowercased, "" if none)
//
// Display inputs (original values — written into output columns as-is):
//   orig_name1/2        — original Name (not lowercased)
//   orig_syn1/2         — original synonyms joined by "; "
//   orig_ik1/2          — original InChIKey
//
// Returns DataFrame with columns:
//   Formula, InChIKey_lib1, InChIKey_lib2,
//   Name_lib1, Synonyms_lib1, Name_lib2, Synonyms_lib2, Match_type
// -----------------------------------------------------------------------------
// [[Rcpp::export]]
DataFrame compare_libraries_cpp(
    CharacterVector formula1,   CharacterVector formula2,
    CharacterVector ik1,        CharacterVector ik2,
    CharacterVector name1,      CharacterVector name2,
    CharacterVector syn1,       CharacterVector syn2,
    CharacterVector orig_name1, CharacterVector orig_name2,
    CharacterVector orig_syn1,  CharacterVector orig_syn2,
    CharacterVector orig_ik1,   CharacterVector orig_ik2
) {
  int n1 = formula1.size();
  int n2 = formula2.size();

  // Build hash map: formula string -> list of lib2 indices with that formula
  std::unordered_map<std::string, std::vector<int>> fmap;
  fmap.reserve(n2 * 2);
  for (int j = 0; j < n2; ++j) {
    std::string f = Rcpp::as<std::string>(formula2[j]);
    if (!f.empty()) fmap[f].push_back(j);
  }

  std::vector<std::string> out_formula, out_ik1, out_ik2;
  std::vector<std::string> out_name1, out_name2;
  std::vector<std::string> out_syn1,  out_syn2;
  std::vector<std::string> out_match;

  for (int i = 0; i < n1; ++i) {
    std::string f1 = Rcpp::as<std::string>(formula1[i]);
    if (f1.empty()) continue;

    auto it = fmap.find(f1);
    if (it == fmap.end()) continue;

    // Cache lib1 strings for this compound
    std::string s_ik1   = Rcpp::as<std::string>(ik1[i]);
    std::string s_name1 = Rcpp::as<std::string>(name1[i]);
    std::string s_syn1  = Rcpp::as<std::string>(syn1[i]);
    std::string d_ik1   = Rcpp::as<std::string>(orig_ik1[i]);
    std::string d_name1 = Rcpp::as<std::string>(orig_name1[i]);
    std::string d_syn1  = Rcpp::as<std::string>(orig_syn1[i]);

    for (int j : it->second) {
      std::string match = classify_pair(
        s_ik1,   Rcpp::as<std::string>(ik2[j]),
        s_name1, Rcpp::as<std::string>(name2[j]),
        s_syn1,  Rcpp::as<std::string>(syn2[j])
      );
      out_formula.push_back(f1);
      out_ik1.push_back(d_ik1);
      out_ik2.push_back(Rcpp::as<std::string>(orig_ik2[j]));
      out_name1.push_back(d_name1);
      out_name2.push_back(Rcpp::as<std::string>(orig_name2[j]));
      out_syn1.push_back(d_syn1);
      out_syn2.push_back(Rcpp::as<std::string>(orig_syn2[j]));
      out_match.push_back(match);
    }
  }

  return DataFrame::create(
    _["Formula"]          = out_formula,
    _["InChIKey_lib1"]    = out_ik1,
    _["InChIKey_lib2"]    = out_ik2,
    _["Name_lib1"]        = out_name1,
    _["Synonyms_lib1"]    = out_syn1,
    _["Name_lib2"]        = out_name2,
    _["Synonyms_lib2"]    = out_syn2,
    _["Match_type"]       = out_match,
    _["stringsAsFactors"] = false
  );
}
