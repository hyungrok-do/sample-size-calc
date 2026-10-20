// Name-based dispatcher for the distribution functions (used by ssc_eval()
// and the reference-value tests).
#include "dispatch.hpp"

#include <unordered_map>

#include "stats.hpp"

namespace ssc::stats {

namespace {
using Fn = double (*)(const double*);
struct Entry {
  int nargs;
  Fn fn;
};

inline bool B(double v) { return v != 0; }

const std::unordered_map<std::string, Entry>& table() {
  static const std::unordered_map<std::string, Entry> t = {
      {"lgamma", {1, [](const double* a) { return lgammafn(a[0]); }}},
      {"lbeta", {2, [](const double* a) { return lbeta(a[0], a[1]); }}},
      {"stirlerr", {1, [](const double* a) { return stirlerr(a[0]); }}},
      {"bd0", {2, [](const double* a) { return bd0(a[0], a[1]); }}},
      {"pnorm", {2, [](const double* a) { return pnorm(a[0], B(a[1])); }}},
      {"dnorm", {1, [](const double* a) { return dnorm(a[0]); }}},
      {"qnorm", {2, [](const double* a) { return qnorm(a[0], B(a[1])); }}},
      {"pgamma", {3, [](const double* a) { return pgamma_raw(a[0], a[1], B(a[2])); }}},
      {"pcorr", {4, [](const double* a) { return pcorr(a[0], a[1], a[2], B(a[3])); }}},
      {"dcorr", {3, [](const double* a) { return dcorr(a[0], a[1], a[2]); }}},
      {"qgamma", {3, [](const double* a) { return qgamma(a[0], a[1], B(a[2])); }}},
      {"pchisq", {3, [](const double* a) { return pchisq(a[0], a[1], B(a[2])); }}},
      {"qchisq", {3, [](const double* a) { return qchisq(a[0], a[1], B(a[2])); }}},
      {"dchisq", {2, [](const double* a) { return dchisq(a[0], a[1]); }}},
      {"dpois", {2, [](const double* a) { return dpois(a[0], a[1]); }}},
      {"ppois", {3, [](const double* a) { return ppois(a[0], a[1], B(a[2])); }}},
      {"qpois", {3, [](const double* a) { return qpois(a[0], a[1], B(a[2])); }}},
      {"pbeta", {4, [](const double* a) { return pbeta(a[0], a[1], a[2], B(a[3])); }}},
      {"dbeta", {3, [](const double* a) { return dbeta(a[0], a[1], a[2], 1 - a[0]); }}},
      {"qbeta", {4, [](const double* a) { return qbeta(a[0], a[1], a[2], B(a[3])); }}},
      {"pt", {3, [](const double* a) { return pt(a[0], a[1], B(a[2])); }}},
      {"dt", {2, [](const double* a) { return dt(a[0], a[1]); }}},
      {"qt", {3, [](const double* a) { return qt(a[0], a[1], B(a[2])); }}},
      {"pf", {4, [](const double* a) { return pf(a[0], a[1], a[2], B(a[3])); }}},
      {"qf", {4, [](const double* a) { return qf(a[0], a[1], a[2], B(a[3])); }}},
      {"dbinom", {3, [](const double* a) { return dbinom(a[0], a[1], a[2]); }}},
      {"pbinom", {4, [](const double* a) { return pbinom(a[0], a[1], a[2], B(a[3])); }}},
      {"qbinom", {3, [](const double* a) { return qbinom(a[0], a[1], a[2]); }}},
      {"dhyper", {4, [](const double* a) { return dhyper(a[0], a[1], a[2], a[3]); }}},
      {"phyper", {5, [](const double* a) { return phyper(a[0], a[1], a[2], a[3], B(a[4])); }}},
      {"pnt", {4, [](const double* a) { return pnt(a[0], a[1], a[2], B(a[3])); }}},
      {"pnchisq", {4, [](const double* a) { return pnchisq(a[0], a[1], a[2], B(a[3])); }}},
      {"pnf", {5, [](const double* a) { return pnf(a[0], a[1], a[2], a[3], B(a[4])); }}},
      {"pnbeta", {5, [](const double* a) { return pnbeta(a[0], a[1], a[2], a[3], B(a[4]), 1 - a[0]); }}},
      {"pbvnorm", {3, [](const double* a) { return pbvnorm(a[0], a[1], a[2]); }}},
  };
  return t;
}
}  // namespace

bool eval(const std::string& name, const double* args, int nargs, double& out) {
  const auto& t = table();
  auto it = t.find(name);
  if (it == t.end() || it->second.nargs != nargs) return false;
  out = it->second.fn(args);
  return true;
}

std::vector<std::string> function_names() {
  std::vector<std::string> v;
  for (const auto& kv : table()) v.push_back(kv.first);
  return v;
}

}  // namespace ssc::stats
