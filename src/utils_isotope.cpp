// utils_isotope_fixed.cpp


#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <string>
#include <numeric>
#include <algorithm>
#include <set>

// ── Inline sort helpers (from utils_sort.h) ──────────────────────
static std::vector<int> compute_sorted_index(const Rcpp::NumericVector& v) {
  std::vector<int> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](int a, int b) {
    return v[a] < v[b];
  });
  return idx;
}

static Rcpp::NumericVector reorder_num(const Rcpp::NumericVector& v,
                                       const std::vector<int>& idx) {
  Rcpp::NumericVector out(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) out[i] = v[idx[i]];
  return out;
}

static Rcpp::CharacterVector reorder_chr(const Rcpp::CharacterVector& v,
                                         const std::vector<int>& idx) {
  Rcpp::CharacterVector out(idx.size());
  for (size_t i = 0; i < idx.size(); ++i) out[i] = v[idx[i]];
  return out;
}
// ─────────────────────────────────────────────────────────────────

// [[Rcpp::export]]
Rcpp::CharacterVector Check_Iso_Charge(
    Rcpp::CharacterVector name,
    Rcpp::NumericVector   mz,
    Rcpp::NumericVector   rt,
    Rcpp::NumericVector   intensity,
    double iso_threshold,
    int    max_charge,
    int    max_isotope,
    int&   groupIndex,
    double ppm          = 5.0,
    double mzabs        = 0.001,
    double rt_threshold = 3.0)
{
  const int n = name.size();

  // ── Sort all inputs by mz ascending ────────────────────────────
  // If FT is RT-ordered; without this sort, the break on
  // mz_diff > window fires on the first high-mz unrelated feature,
  // skipping any isotope partner that appears later in the array.
  std::vector<int> idx = compute_sorted_index(mz);
  Rcpp::CharacterVector name_s = reorder_chr(name,      idx);
  Rcpp::NumericVector   mz_s   = reorder_num(mz,        idx);
  Rcpp::NumericVector   rt_s   = reorder_num(rt,        idx);
  Rcpp::NumericVector   int_s  = reorder_num(intensity, idx);
  // ───────────────────────────────────────────────────────────────

  Rcpp::CharacterVector labels(n, "");   // in sorted space

  for (int i = 0; i < n; ++i) {

    if (labels[i] != "") continue;

    bool foundAny   = false;
    int  confirmedC = -1;

    for (int c = 1; c <= max_charge; ++c) {

      if (confirmedC != -1) break;   // lowest confirmed charge wins

      int maxIsoFound = 0;

      // claimed_m: tracks which isotope positions (m=1,2,...) have already
      // been assigned a label in this c-pass.  Because the j-loop walks in
      // mz-sorted ascending order, the first j that passes all checks for a
      // given m is always the nearest-mz match — exactly what we want.
      // Any subsequent j that also matches the same m (e.g. FT0231 at +10 mDa
      // after FT0230 already claimed [M+2]+) is skipped for the label write
      // but still advances maxIsoFound so higher isotopes remain reachable.
      std::set<int> claimed_m;

      for (int j = i + 1; j < n; ++j) {

        // NOTE: do NOT skip already-labeled j here.
        // If we skip the entire j when labels[j] != "", maxIsoFound never
        // advances past the already-labeled position, and the sequential
        // enforcement (m > maxIsoFound + 1) then gates out all higher
        // isotopes even though a valid m=1 partner exists at that position.
        // Instead we let the m-loop run normally and only block the label
        // write below (not the maxIsoFound advance).

        double mz_diff = mz_s[j] - mz_s[i];   // always >= 0 after sort
        double rt_diff = std::fabs(rt_s[j] - rt_s[i]);

        // Early exit: mz_diff already past the largest possible isotope peak
        if (mz_diff > (max_isotope * iso_threshold / c + 0.02)) break;

        for (int m = 1; m <= max_isotope; ++m) {

          if (m > maxIsoFound + 1) continue;   // sequential enforcement

          double expected = m * iso_threshold / c;
          double tol = (m == 1)
            ? (mz_s[i] * ppm * 1e-6 + mzabs) / c
          : 0.010;

          if (std::fabs(mz_diff - expected) > tol) continue;
          if (rt_diff >= rt_threshold) continue;

          // Intensity checks
          //
          // m=1, c=1: [M+1]+ is the C13 isotope peak.
          //   Expected ratio = (neutral_mass / 12) * 1.1%
          //   For singly charged: neutral_mass ≈ mz - 1.00728 (H)
          //   We apply a WINDOW check (not just upper bound):
          //     min = expected * 0.05  (allow very low — detection noise)
          //     max = expected * 4.0   (allow up to 4x — adduct overlap,
          //                             co-eluting isobaric species)
          //   This correctly handles large lipids (mz~700+) where a naive
          //   upper-bound of (mz/12)*0.011 ≈ 64% is far too permissive
          //   and would accept unrelated co-eluting features as satellites.
          //   A window forces the satellite to be in a physically realistic
          //   range around the expected C13 abundance.
          //
          // m=1, c>=2: [M+1]2+ etc. — isotope envelope shape depends on
          //   charge and mass.  Apply same window logic with expected ratio
          //   computed from neutral mass = (mz * c) - c * 1.00728.
          //
          // m>=2, any c: must be <= [M+1] intensity (isotope envelope
          //   decreases or plateaus for common metabolomics elements).
          //   Use loose upper bound of monoisotope intensity.
          if (m == 1) {
            // Neutral mass approximation: mz * c - c * 1.00728
            double neutral_mass = mz_s[i] * c - c * 1.00728;
            if (neutral_mass < 1.0) neutral_mass = mz_s[i];  // fallback
            double expected_ratio = (neutral_mass / 12.0) * 0.011;
            double min_ratio = expected_ratio * 0.05;
            double max_ratio = expected_ratio * 4.0;
            double obs_ratio = (int_s[i] > 0.0)
              ? int_s[j] / int_s[i] : 0.0;
            if (obs_ratio < min_ratio || obs_ratio > max_ratio) continue;
          } else {
            // m>=2: satellite must not exceed monoisotope
            if (int_s[j] > int_s[i]) continue;
          }

          // ── All checks passed ───────────────────────────────
          // Label monoisotope i on first confirmed satellite
          if (confirmedC == -1) {
            std::string lbl = "[" + std::to_string(groupIndex) + "][M]";
            if (c > 1) lbl += std::to_string(c);
            lbl += "+";
            labels[i]  = Rcpp::String(lbl);
            confirmedC = c;
            foundAny   = true;
          }

          // Write satellite label only if:
          //   (a) j not yet labeled by any earlier monoisotope, AND
          //   (b) this isotope position m not yet claimed in this c-pass.
          // Condition (b) ensures nearest-mz wins: the first j that reaches
          // this point for a given m is always the closest one (mz-sorted).
          // ALWAYS advance maxIsoFound regardless, so m+1 remains reachable.
          if (labels[j] == "" && claimed_m.find(m) == claimed_m.end()) {
            std::string lbl = "[" + std::to_string(groupIndex) +
              "][M+" + std::to_string(m) + "]";
            if (c > 1) lbl += std::to_string(c);
            lbl += "+";
            labels[j] = Rcpp::String(lbl);
            claimed_m.insert(m);
          }

          maxIsoFound = m;   // advance regardless of label write outcome

        } // m
      }   // j
    }     // c

    if (foundAny) groupIndex++;

  } // i

  // ── Unsort: map sorted labels back to original input order ─────
  // idx[sorted_pos] = original_pos, so writing labels_out[idx[i]]
  // inverts the permutation correctly.
  Rcpp::CharacterVector labels_out(n, "");
  for (int i = 0; i < n; ++i) {
    labels_out[idx[i]] = labels[i];
  }
  return labels_out;
}



// #include <Rcpp.h>
// #include <cmath>
// #include <vector>
// #include <string>
// #include <numeric>
// #include <algorithm>
//
// // ── Inline sort helpers (from utils_sort.h) ──────────────────────
// static std::vector<int> compute_sorted_index(const Rcpp::NumericVector& v) {
//   std::vector<int> idx(v.size());
//   std::iota(idx.begin(), idx.end(), 0);
//   std::sort(idx.begin(), idx.end(), [&](int a, int b) {
//     return v[a] < v[b];
//   });
//   return idx;
// }
//
// static Rcpp::NumericVector reorder_num(const Rcpp::NumericVector& v,
//                                        const std::vector<int>& idx) {
//   Rcpp::NumericVector out(idx.size());
//   for (size_t i = 0; i < idx.size(); ++i) out[i] = v[idx[i]];
//   return out;
// }
//
// static Rcpp::CharacterVector reorder_chr(const Rcpp::CharacterVector& v,
//                                          const std::vector<int>& idx) {
//   Rcpp::CharacterVector out(idx.size());
//   for (size_t i = 0; i < idx.size(); ++i) out[i] = v[idx[i]];
//   return out;
// }
// // ─────────────────────────────────────────────────────────────────
//
// // [[Rcpp::export]]
// Rcpp::CharacterVector Check_Iso_Charge(
//     Rcpp::CharacterVector name,
//     Rcpp::NumericVector   mz,
//     Rcpp::NumericVector   rt,
//     Rcpp::NumericVector   intensity,
//     double iso_threshold,
//     int    max_charge,
//     int    max_isotope,
//     int&   groupIndex,
//     double ppm   = 5.0,
//     double mzabs = 0.001)
// {
//   const int n = name.size();
//
//   // ── Sort all inputs by mz ascending ────────────────────────────
//   // xcms output is RT-ordered; without this sort, the break on
//   // mz_diff > window fires on the first high-mz unrelated feature,
//   // skipping any isotope partner that appears later in the array.
//   std::vector<int> idx = compute_sorted_index(mz);
//   Rcpp::CharacterVector name_s = reorder_chr(name,      idx);
//   Rcpp::NumericVector   mz_s   = reorder_num(mz,        idx);
//   Rcpp::NumericVector   rt_s   = reorder_num(rt,        idx);
//   Rcpp::NumericVector   int_s  = reorder_num(intensity, idx);
//   // ───────────────────────────────────────────────────────────────
//
//   Rcpp::CharacterVector labels(n, "");   // in sorted space
//
//   for (int i = 0; i < n; ++i) {
//
//     if (labels[i] != "") continue;
//
//     bool foundAny   = false;
//     int  confirmedC = -1;
//
//     for (int c = 1; c <= max_charge; ++c) {
//
//       if (confirmedC != -1) break;   // lowest confirmed charge wins
//
//       int maxIsoFound = 0;
//
//       for (int j = i + 1; j < n; ++j) {
//
//         if (labels[j] != "") continue;
//
//         double mz_diff = mz_s[j] - mz_s[i];   // always >= 0 after sort
//         double rt_diff = std::fabs(rt_s[j] - rt_s[i]);
//
//         // Early exit: mz_diff already past the largest possible isotope peak
//         if (mz_diff > (max_isotope * iso_threshold / c + 0.02)) break;
//
//         for (int m = 1; m <= max_isotope; ++m) {
//
//           if (m > maxIsoFound + 1) continue;   // sequential enforcement
//
//           double expected = m * iso_threshold / c;
//           double tol = (m == 1)
//             ? (mz_s[i] * ppm * 1e-6 + mzabs) / c
//           : 0.010;
//
//           if (std::fabs(mz_diff - expected) > tol) continue;
//           if (rt_diff >= 3.0) continue;
//
//           // Intensity checks
//           //
//           // m=1, c=1: [M+1]+ is the C13 peak. Intensity is
//           //   predictable: expected = intensity[M] * (mass/12) * 0.011
//           //   Use as strict upper bound (C13 abundance well-defined).
//           //   Note: mz_s[i] ≈ mass for singly charged, so no *c needed.
//           //
//           // m=1, c>=2: [M+1]c+ sits at +1/c Da in m/z space, i.e.
//           //   the fractional-Da isotope peak of the multiply charged
//           //   envelope. Its intensity is NOT governed by C13 alone —
//           //   it depends on the full isotope envelope shape sampled at
//           //   1/c Da spacing. Can legitimately be 30-60% of [M]c+
//           //   (e.g. charge=2 at mz~90 -> mass~178, [M+1]/[M] ~40%).
//           //   Only require it does not exceed monoisotopic.
//           //
//           // m>=2, any c: require <= monoisotopic (conservative but
//           //   correct for all common metabolomics elements).
//           if (m == 1 && c == 1) {
//             double max_ratio = (mz_s[i] / 12.0) * 0.011;
//             if (int_s[j] > int_s[i] * max_ratio) continue;
//           } else {
//             if (int_s[j] > int_s[i]) continue;
//           }
//
//           // ── All checks passed ───────────────────────────────
//           if (confirmedC == -1) {
//             std::string lbl = "[" + std::to_string(groupIndex) + "][M]";
//             if (c > 1) lbl += std::to_string(c);
//             lbl += "+";
//             labels[i]  = Rcpp::String(lbl);
//             confirmedC = c;
//             foundAny   = true;
//           }
//
//           if (labels[j] == "") {
//             std::string lbl = "[" + std::to_string(groupIndex) +
//               "][M+" + std::to_string(m) + "]";
//             if (c > 1) lbl += std::to_string(c);
//             lbl += "+";
//             labels[j] = Rcpp::String(lbl);
//           }
//
//           maxIsoFound = m;
//
//         } // m
//       }   // j
//     }     // c
//
//     if (foundAny) groupIndex++;
//
//   } // i
//
//   // ── Unsort: map sorted labels back to original input order ─────
//   // idx[sorted_pos] = original_pos, so writing labels_out[idx[i]]
//   // inverts the permutation correctly.
//   Rcpp::CharacterVector labels_out(n, "");
//   for (int i = 0; i < n; ++i) {
//     labels_out[idx[i]] = labels[i];
//   }
//   return labels_out;
// }
