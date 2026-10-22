// Name-based dispatcher for the distribution functions (used by ssc_eval()
// and the reference-value tests).
#include "dispatch.hpp"

#include <cmath>
#include <limits>
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

// Whether the arguments lie in the function's domain, as R checks them: the internal routines assume it (a
// negative shape or an infinite count can make them loop for ever), and the calculations always respect it.
bool in_domain(const std::string& name, const double* a, int nargs) {
  for (int i = 0; i < nargs; ++i)
    if (std::isnan(a[i])) return false;
  auto pos = [](double v) { return v > 0 && std::isfinite(v); };       // a shape or degrees of freedom
  auto count = [](double v) { return v >= 0 && std::isfinite(v); };    // a count, a rate or a noncentrality
  auto prob = [](double v) { return v >= 0 && v <= 1; };
  auto corr = [](double v) { return v >= -1 && v <= 1; };
  if (name == "stirlerr") return count(a[0]);
  if (name == "bd0") return count(a[0]) && count(a[1]);
  if (name == "pgamma") return pos(a[1]);
  if (name == "qgamma") return prob(a[0]) && pos(a[1]);
  if (name == "pchisq" || name == "dchisq") return pos(a[1]);
  if (name == "qchisq") return prob(a[0]) && pos(a[1]);
  if (name == "dpois" || name == "ppois") return count(a[1]);
  if (name == "qpois") return prob(a[0]) && count(a[1]);
  if (name == "pbeta" || name == "dbeta") return pos(a[1]) && pos(a[2]);
  if (name == "qbeta") return prob(a[0]) && pos(a[1]) && pos(a[2]);
  if (name == "pt" || name == "dt") return a[1] > 0;  // infinite degrees of freedom: the normal distribution
  if (name == "qt") return prob(a[0]) && a[1] > 0;
  if (name == "pf") return a[1] > 0 && a[2] > 0;
  if (name == "qf") return prob(a[0]) && a[1] > 0 && a[2] > 0;
  if (name == "dbinom" || name == "pbinom") return count(a[1]) && prob(a[2]);
  if (name == "qbinom") return prob(a[0]) && count(a[1]) && prob(a[2]);
  if (name == "dhyper" || name == "phyper") return count(a[1]) && count(a[2]) && count(a[3]);
  if (name == "pnt") return a[1] > 0 && std::isfinite(a[2]);
  if (name == "pnchisq") return pos(a[1]) && count(a[2]);
  if (name == "pnf") return a[1] > 0 && a[2] > 0 && count(a[3]);
  if (name == "pnbeta") return pos(a[1]) && pos(a[2]) && count(a[3]);
  if (name == "pcorr" || name == "dcorr") return corr(a[1]) && pos(a[2]);
  if (name == "pbvnorm") return corr(a[2]);
  if (name == "qnorm") return prob(a[0]);
  return true;
}
}  // namespace

bool eval(const std::string& name, const double* args, int nargs, double& out) {
  const auto& t = table();
  auto it = t.find(name);
  if (it == t.end() || it->second.nargs != nargs) return false;
  out = in_domain(name, args, nargs) ? it->second.fn(args) : std::numeric_limits<double>::quiet_NaN();
  return true;
}

std::vector<std::string> function_names() {
  std::vector<std::string> v;
  for (const auto& kv : table()) v.push_back(kv.first);
  return v;
}

}  // namespace ssc::stats
