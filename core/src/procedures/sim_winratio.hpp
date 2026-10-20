// The pairwise comparisons of the win-ratio simulation (sim_winratio.cpp), which the tests check directly:
// a fast exact count by sorting and a direct count over every pair.
#pragma once

#include <vector>

namespace ssc::procs::winratio {

/// One patient as observed: the death time if it happened during follow-up (else ∞) and the time until which
/// the patient is known to be alive (death or end of follow-up); the first non-fatal event if it happened
/// during follow-up and before death (else ∞) and the time until which the patient is known to be free of it;
/// the measured outcome if alive at the end of follow-up (else NaN).  A priority not in use has ∞ times (or
/// NaN values) for everyone.
struct Patient {
  double death, alive, event, free, value;
};

/// +1 if patient a wins against b, −1 if a loses, 0 for a tie, over the time both were followed: a loses on
/// death when a died while b was still known to be alive (Gehan's rule), and likewise for the non-fatal event;
/// a measured outcome decides only between two patients alive at the end of their follow-up, by more than the
/// margin.
inline int compare(const Patient& a, const Patient& b, double margin) {
  const int s1 = (b.death < a.alive) - (a.death < b.alive);
  const int s2 = (b.event < a.free) - (a.event < b.free);
  const int s3 = (b.value < a.value - margin) - (b.value > a.value + margin);  // NaN compares false: a tie
  return s1 != 0 ? s1 : (s2 != 0 ? s2 : s3);
}

/// For every patient k, the number of patients of group g (0: the first n1 patients, 1: the others) that k
/// beats (win[g][k]) and that beat k (loss[g][k]).
struct Counts {
  std::vector<double> win[2], loss[2];
};

/// Exact counts in O(N log N): each priority's condition is a dominance count, done by sweeping sorted
/// patients through a Fenwick tree.
void pair_counts(const std::vector<Patient>& patients, long n1, double margin, Counts& out);

/// The same counts by comparing every pair (O(N²)), for checking pair_counts.
void pair_counts_direct(const std::vector<Patient>& patients, long n1, double margin, Counts& out);

}  // namespace ssc::procs::winratio
