#include "model.hpp"
#include "study.hpp"
#include "sim/random.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>

namespace ssc::dag {
[[noreturn]] void fail(const std::string &message) {
  throw Error("invalid_design", message);
}
double number(const J &o, const std::string &key) {
  const J &v = o.get(key);
  if (!v.is_number() || !std::isfinite(v.as_number()))
    fail(key + " must be a finite number.");
  return v.as_number();
}
int integer(const J &o, const std::string &key, int lo, int hi, int fallback) {
  if (!o.find(key))
    return fallback;
  const double v = number(o, key);
  if (v != std::floor(v) || v < lo || v > hi)
    fail(key + " is outside its supported integer range.");
  return static_cast<int>(v);
}
J object(std::initializer_list<std::pair<std::string, J>> values) {
  J out = J::object();
  for (const auto &kv : values)
    out.set(kv.first, kv.second);
  return out;
}
J array(const Vec &values) {
  J out = J::array();
  for (double v : values)
    out.push(v);
  return out;
}
int Graph::index(const std::string &id) const {
  for (std::size_t i = 0; i < nodes.size(); ++i)
    if (nodes[i].id == id)
      return static_cast<int>(i);
  fail("An arrow, intervention or estimand requires an existing variable: " +
       id);
}
static bool identifier(const std::string &s) {
  const auto letter = [](char c) {
    return (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z');
  };
  if (s.empty() || s.size() > 32 || !letter(s[0]))
    return false;
  for (char c : s)
    if (!letter(c) && !(c >= '0' && c <= '9') && c != '_')
      return false;
  return true;
}
Graph parse_graph(const J &raw) {
  if (raw.get("version").as_number() != 1)
    fail("Unsupported DAG file version.");
  const auto &nodes = raw.get("nodes");
  const auto &edges = raw.get("edges");
  if (!nodes.is_array() || nodes.size() > 40 || !edges.is_array() ||
      edges.size() > 780)
    fail("Invalid DAG size (at most 40 variables).");
  Graph g;
  g.json = raw;
  std::set<std::string> ids;
  const std::map<std::string, std::string> kinds = {
      {"gaussian", "continuous"},  {"uniform", "continuous"},
      {"lognormal", "continuous"}, {"gamma", "continuous"},
      {"student", "continuous"},   {"beta", "continuous"},
      {"bernoulli", "binary"},     {"poisson", "count"},
      {"negbin", "count"},
      {"weibull", "survival"},     {"exponential", "survival"}};
  for (const auto &v : nodes.elements()) {
    Node node;
    node.id = v.get("id").as_string();
    node.label = v.get("label").as_string();
    if (!identifier(node.id) || !ids.insert(node.id).second)
      fail("Variable IDs must be unique valid identifiers.");
    if (node.label.empty() || node.label.size() > 320)
      fail("Each variable needs a label.");
    node.kind = v.get("kind").as_string();
    node.family = v.get("family").as_string();
    if (!kinds.count(node.family) || kinds.at(node.family) != node.kind)
      fail("Choose a distribution compatible with the variable type.");
    node.intercept = number(v, "intercept");
    auto positive = [&](const char *key) {
      const double x = number(v, key);
      if (!(x > 0))
        fail(std::string(key) + " must be positive.");
      return x;
    };
    if (node.family == "gaussian" || node.family == "lognormal" || node.family == "student")
      node.sd = positive("sd");
    if (node.family == "student" && !((node.df = positive("df")) > 1))
      fail("A Student t variable needs more than 1 degree of freedom, so that its mean exists.");
    if (node.family == "beta")
      node.precision = positive("precision");
    if (node.family == "negbin")
      node.overdispersion = positive("overdispersion");
    if (v.find("measured")) {
      if (!v.get("measured").is_bool())
        fail("measured must be true or false.");
      node.measured = v.get("measured").as_bool();
    }
    if (node.family == "uniform")
      node.width = positive("width");
    if (node.family == "gamma" || node.family == "weibull")
      node.shape = positive("shape");
    if (node.kind == "survival") {
      node.scale = positive("scale");
      const auto &c = v.get("censoring");
      node.censor = c.get("mode").as_string();
      if (node.censor != "none" && node.censor != "administrative" &&
          node.censor != "variable")
        fail("Choose a censoring observation rule.");
      if (node.censor == "administrative") {
        node.cutoff = number(c, "time");
        if (!(node.cutoff > 0))
          fail("The administrative cutoff must be positive.");
      }
    }
    // Stable identifier hash, never std::hash (implementation dependent).
    node.stream = 14695981039346656037ULL;
    for (unsigned char c : node.id) {
      node.stream ^= c;
      node.stream *= 1099511628211ULL;
    }
    g.nodes.push_back(std::move(node));
  }
  std::vector<int> indegree(g.nodes.size());
  std::vector<std::vector<int>> children(g.nodes.size());
  std::set<std::pair<int, int>> pairs;
  for (const auto &e : edges.elements()) {
    const int source = g.index(e.get("source").as_string()),
              target = g.index(e.get("target").as_string());
    if (source == target || !pairs.insert({source, target}).second)
      fail("Duplicate or self-referencing arrow.");
    g.nodes[target].terms.push_back({source, -1, number(e, "coefficient")});
    ++indegree[target];
    children[source].push_back(target);
  }
  for (int i = 0; i < static_cast<int>(g.nodes.size()); ++i)
    if (!indegree[i])
      g.order.push_back(i);
  for (std::size_t at = 0; at < g.order.size(); ++at)
    for (int child : children[g.order[at]])
      if (--indegree[child] == 0)
        g.order.push_back(child);
  if (g.order.size() != g.nodes.size())
    fail("That arrow creates a cycle. A DAG must have no directed cycles.");
  for (int i = 0; i < static_cast<int>(g.nodes.size()); ++i) {
    Node &node = g.nodes[i];
    const auto &v = nodes.at(i);
    const auto &terms = v.get("interactions");
    if (!terms.is_array() || terms.size() > 1600)
      fail("Invalid interaction list.");
    std::set<std::pair<int, int>> interactions;
    for (const auto &t : terms.elements()) {
      const int left = g.index(t.get("left").as_string()),
                right = g.index(t.get("right").as_string());
      if (left == right || !pairs.count({left, i}) || !pairs.count({right, i}))
        fail("Interactions require two different parents.");
      if (!interactions.insert(std::minmax(left, right)).second)
        fail("That interaction already exists.");
      node.terms.push_back({left, right, number(t, "coefficient")});
    }
    if (v.find("functions")) {
      if (!v.get("functions").is_array() || v.get("functions").size() > 120)
        fail("Use at most 120 additional response terms per variable.");
      std::set<std::string> seen;
      for(const auto &t:node.terms) seen.insert(term_key(t));
      for (const auto &value : v.get("functions").elements()) {
        const auto t = parse_term(g, value);
        if (!pairs.count({t.left, i}) || (t.right >= 0 && !pairs.count({t.right, i})))
          fail("Response terms must use existing parents of the child.");
        if (!seen.insert(term_key(t)).second) fail("Duplicate response term.");
        node.terms.push_back(t);
      }
    }
    if (node.kind != "survival")
      continue;
    if (node.censor == "variable") {
      node.censor_index =
          g.index(v.get("censoring").get("variable").as_string());
      if (node.censor_index == i)
        fail("Select a different variable for censoring.");
    }
    for (const auto &suffix : {"_time", "_event"})
      if (!ids.insert(node.id + suffix).second)
        fail("A derived survival column conflicts with a variable ID.");
  }
  if(raw.get("study").is_null()) for(const auto &v:nodes.elements()) {
    if((v.find("scope") && v.get("scope").as_string()!="subject") || v.get("lags").size()) fail("Enable study structure before setting repeated or cluster variables.");
    for(const auto &kv:v.get("repeated").items()) if(number(v.get("repeated"),kv.first)!=0) fail("Repeated effects require study structure.");
  }
  parse_study(g);
  return g;
}
static double poisson(double mean, sim::Rng &rng) {
  // Exact Poisson additivity, bounded means; no normal approximation.
  double total = 0;
  for (double remaining = mean; remaining > 0;) {
    const double chunk = std::min(remaining, 20.0), limit = std::exp(-chunk);
    double product = 1;
    int k = 0;
    do {
      ++k;
      product *= rng.uniform();
    } while (product > limit);
    total += k - 1;
    remaining -= chunk;
  }
  return total;
}
double draw_value(const Node &node, double eta, sim::Rng &rng, double normal) {
  if(!std::isfinite(eta)) throw Error("generation_error", "A conditional predictor overflowed.");
  double value;
      if (node.family == "gaussian")
        value = eta + node.sd * (std::isnan(normal) ? rng.normal() : normal);
      else if (node.family == "uniform")
        value = eta + node.width * (rng.uniform() - .5);
      else if (node.family == "lognormal")
        value = std::exp(eta + node.sd * (std::isnan(normal) ? rng.normal() : normal));
      else if (node.family == "gamma")
        value = rng.gamma(node.shape) / node.shape * std::exp(eta);
      else if (node.family == "student")
        value = eta + node.sd * (std::isnan(normal) ? rng.normal() : normal) /
                          std::sqrt(2 * rng.gamma(node.df / 2) / node.df);
      else if (node.family == "beta") {
        // Mean expit(eta), variance mean (1 - mean) / (1 + precision).
        const double mean = eta >= 0 ? 1 / (1 + std::exp(-eta)) : std::exp(eta) / (1 + std::exp(eta));
        const double a = rng.gamma(mean * node.precision), b = rng.gamma((1 - mean) * node.precision);
        // Shapes far below 1 can underflow a gamma draw to zero; keep the
        // value strictly inside (0, 1).
        value = std::clamp(a + b > 0 ? a / (a + b) : mean, 1e-300, 1 - 0x1p-53);
      }
      else if (node.family == "bernoulli")
        value = rng.uniform() < (eta >= 0 ? 1 / (1 + std::exp(-eta))
                                          : std::exp(eta) / (1 + std::exp(eta)))
                    ? 1
                    : 0;
      else if (node.family == "poisson") {
        const double mean = std::exp(eta);
        if (!std::isfinite(mean) || !(mean > 0) || mean > 1000)
          throw Error("generation_error",
                      node.label + ": the Poisson rate must be in (0, 1000].");
        value = poisson(mean, rng);
      } else if (node.family == "negbin") {
        // Gamma-Poisson mixture: a rate with mean exp(eta) and variance
        // overdispersion * exp(2 eta).
        const double mean = std::exp(eta), k = 1 / node.overdispersion;
        if (!std::isfinite(mean) || !(mean > 0) || mean > 1000)
          throw Error("generation_error",
                      node.label + ": the negative binomial mean must be in (0, 1000].");
        const double rate = mean * rng.gamma(k) / k;
        if (!(rate <= 100000))
          throw Error("generation_error",
                      node.label + ": a negative binomial rate exceeded 100,000; reduce the overdispersion.");
        value = rate > 0 ? poisson(rate, rng) : 0;
      } else
        value = std::exp(std::log(node.scale) +
                         (std::log(rng.exponential()) - eta) / node.shape);
      if (!std::isfinite(value) ||
          ((node.family == "lognormal" || node.family == "gamma" ||
            node.family == "beta" || node.kind == "survival") &&
           value <= 0))
        throw Error("generation_error",
                    node.label + ": a simulated value is outside the supported "
                                 "numerical range.");
  return value;
}
Data simulate(const Graph &g, int n, std::uint64_t seed,
              const std::map<int, double> &interventions, bool observe,
              const Token *token) {
  if (g.study.enabled) return simulate_study(g,n,seed,interventions,observe,token);
  if (g.nodes.empty())
    fail("Add at least one variable before simulating data.");
  if (n < 1 || n > 50000)
    fail("Sample size must be a whole number from 1 to 50,000.");
  Data data;
  data.n = n;
  data.values.resize(g.nodes.size(), Vec(n));
  data.times.resize(g.nodes.size());
  data.events.resize(g.nodes.size());
  for (const auto &kv : interventions) {
    const auto &node = g.nodes.at(kv.first);
    const double v = kv.second;
    if (!std::isfinite(v))
      fail("An intervention requires a finite value.");
    if (node.kind == "binary" && v != 0 && v != 1)
      fail("A binary intervention must use 0 or 1.");
    if (node.kind == "count" && (v < 0 || v != std::floor(v)))
      fail("A count intervention must be a nonnegative integer.");
    if (node.kind == "survival" && v <= 0)
      fail("A survival-time intervention must be positive.");
  }
  for (int j : g.order) {
    const auto &node = g.nodes[j];
    auto &column = data.values[j];
    if (interventions.count(j)) {
      std::fill(column.begin(), column.end(), interventions.at(j));
      continue;
    }
    for (int i = 0; i < n; ++i) {
      if (token && (i & 255) == 0)
        token->check();
      // A node and subject have their own calculator RNG stream. Interventions
      // do not shift other draws.
      sim::Rng rng(seed, node.stream, static_cast<std::uint64_t>(i));
      double eta = node.intercept;
      for (const auto &t : node.terms)
        eta += t.coefficient * term_value(t, data.values[t.left][i],
               t.right < 0 ? 1 : data.values[t.right][i]);
      if (!std::isfinite(eta))
        throw Error("generation_error",
                    node.label + ": a conditional predictor overflowed.");
      const double value = draw_value(node, eta, rng);
      column[i] = value;
    }
  }
  if (observe)
    for (int j = 0; j < static_cast<int>(g.nodes.size()); ++j) {
      const auto &node = g.nodes[j];
      if (node.kind != "survival")
        continue;
      data.times[j].resize(n);
      data.events[j].resize(n);
      for (int i = 0; i < n; ++i) {
        const double c = node.censor == "none" ? INFINITY
                         : node.censor == "administrative"
                             ? node.cutoff
                             : data.values[node.censor_index][i];
        if (!(c >= 0))
          throw Error("generation_error",
                      node.label + ": censoring times must be nonnegative.");
        data.times[j][i] = std::min(data.values[j][i], c);
        data.events[j][i] = data.values[j][i] <= c ? 1 : 0;
      }
    }
  return data;
}
J data_json(const Graph &g, const Data &data, std::uint64_t seed,
            bool summarize) {
  J columns = J::object(), summaries = J::array(), survival = J::array(),
    order = J::array();
  auto add = [&](const std::string &id, const Vec &values) {
    columns.set(id, array(values));
    if (!summarize)
      return;
    double mean = 0, m2 = 0, lo = INFINITY, hi = -INFINITY;
    for (int i = 0; i < data.n; ++i) {
      const double v = values[i], delta = v - mean;
      mean += delta / (i + 1);
      m2 += delta * (v - mean);
      lo = std::min(lo, v);
      hi = std::max(hi, v);
    }
    if (!std::isfinite(mean) || !std::isfinite(m2))
      throw Error("generation_error", id + ": summary statistics overflowed.");
    summaries.push(object(
        {{"id", id},
         {"n", data.n},
         {"mean", mean},
         {"sd",
          data.n > 1 ? J(std::sqrt(std::max(0.0, m2 / (data.n - 1)))) : J()},
         {"min", lo},
         {"max", hi}}));
  };
  for (int j = 0; j < static_cast<int>(g.nodes.size()); ++j)
    add(g.nodes[j].id, data.values[j]);
  for (int j = 0; j < static_cast<int>(g.nodes.size()); ++j)
    if (!data.times[j].empty()) {
      add(g.nodes[j].id + "_time", data.times[j]);
      add(g.nodes[j].id + "_event", data.events[j]);
      int count = 0;
      for (double v : data.events[j])
        count += static_cast<int>(v);
      survival.push(object({{"id", g.nodes[j].id},
                            {"events", count},
                            {"censored", data.n - count}}));
    }
  if(g.study.enabled) {
    auto ids=[&](const char *name,const std::vector<int> &x) { Vec v(x.begin(),x.end()); columns.set(name,array(v)); };
    ids("subject_id",data.subject); ids("cluster_id",data.cluster); ids("visit_index",data.visit);
    Vec time; for(int v:data.visit) time.push_back(g.study.times[v]); columns.set("visit_time",array(time));
    for(int j=0;j<static_cast<int>(g.nodes.size());++j) if(g.nodes[j].scope=="visit") {
      columns.set(g.nodes[j].id+"_observed",array(data.observed[j]));
      J recorded=J::array(); for(int r=0;r<data.n;++r) recorded.push(data.observed[j][r] ? J(data.values[j][r]) : J());
      columns.set(g.nodes[j].id+"_recorded",recorded);
    }
  }
  for (int j : g.order)
    order.push(g.nodes[j].id);
  return object({{"n", data.n},
                 {"units",g.study.enabled ? J(data.units) : J(data.n)},
                 {"subjects",g.study.enabled ? J(data.subjects) : J(data.n)},
                 {"sizeUnit",g.study.enabled && !g.study.clusterSizes.empty() ? "clusters" : "subjects"},
                 {"seed", static_cast<unsigned long long>(seed)},
                 {"columns", columns},
                 {"summaries", summaries},
                 {"survival", survival},
                 {"order", order}});
}
Data read_data(const Graph &g, const J &value) {
  Data d;
  d.n = integer(value, "n", 1, 50000, 0);
  if (!d.n)
    fail("A dataset needs its row count.");
  auto column = [&](const std::string &id) {
    const auto &a = value.get("columns").get(id);
    if (!a.is_array() || a.size() != static_cast<std::size_t>(d.n))
      fail("Missing or wrong-length data column: " + id);
    Vec out;
    out.reserve(d.n);
    for (const auto &v : a.elements()) {
      if (!v.is_number() || !std::isfinite(v.as_number()))
        fail("Nonfinite data in " + id);
      out.push_back(v.as_number());
    }
    return out;
  };
  for (const auto &node : g.nodes) {
    d.values.push_back(column(node.id));
    if (node.kind == "survival") {
      d.times.push_back(column(node.id + "_time"));
      d.events.push_back(column(node.id + "_event"));
      for (int i = 0; i < d.n; ++i)
        if (d.times.back()[i] < 0 ||
            (d.events.back()[i] != 0 && d.events.back()[i] != 1))
          fail("Invalid observed survival data.");
    } else {
      d.times.emplace_back();
      d.events.emplace_back();
    }
  }
  if(g.study.enabled) read_study_data(g,value,d);
  return d;
}
std::vector<bool> ancestors(const Graph &g, int node) {
  std::vector<bool> found(g.nodes.size()); found[node]=true; bool changed=true;
  while(changed) { changed=false; for(int j=0;j<static_cast<int>(g.nodes.size());++j) if(found[j]) {
    auto add=[&](int i) { if(i>=0 && !found[i]) { found[i]=true; changed=true; } };
    for(const auto &t:g.nodes[j].terms) { add(t.left); add(t.right); }
    for(const auto &l:g.nodes[j].lags) add(l.source);
  } }
  return found;
}
std::vector<Term> measured_terms(const Graph &g, const std::vector<Term> &terms) {
  std::vector<Term> out;
  std::set<std::string> seen;
  const auto add = [&](const Term &t) {
    if (seen.insert(term_key(t)).second) out.push_back(t);
  };
  // The measured variables an unmeasured one is generated from, nearest first.
  const auto stand_ins = [&](int v) {
    std::vector<int> found, queue{v};
    std::vector<bool> visited(g.nodes.size());
    visited[v] = true;
    for (std::size_t i = 0; i < queue.size(); ++i) {
      std::vector<int> parents;
      for (const auto &t : g.nodes[queue[i]].terms) { parents.push_back(t.left); if (t.right >= 0) parents.push_back(t.right); }
      for (const auto &l : g.nodes[queue[i]].lags) parents.push_back(l.source);
      for (int p : parents) {
        if (visited[p]) continue;
        visited[p] = true;
        if (g.nodes[p].measured) found.push_back(p); else queue.push_back(p);
      }
    }
    std::sort(found.begin(), found.end());
    return found;
  };
  for (const auto &t : terms) {
    const bool left = g.nodes[t.left].measured, right = t.right < 0 || g.nodes[t.right].measured;
    if (left && right) { add(t); continue; }
    for (int v : {t.left, t.right}) {
      if (v < 0) continue;
      if (g.nodes[v].measured) add(Term{v, -1, 0});
      else for (int p : stand_ins(v)) add(Term{p, -1, 0});
    }
  }
  return out;
}
std::vector<bool> descendants(const Graph &g, int node) {
  std::vector<bool> found(g.nodes.size()); found[node]=true; bool changed=true;
  while(changed) { changed=false; for(int j=0;j<static_cast<int>(g.nodes.size());++j) if(!found[j]) {
    bool yes=false;
    for(const auto &t:g.nodes[j].terms) yes |= found[t.left] || (t.right>=0 && found[t.right]);
    for(const auto &l:g.nodes[j].lags) yes |= found[l.source];
    if(yes) { found[j]=true; changed=true; }
  } }
  return found;
}
std::vector<int> calibration_edges(const Graph &g, const J &e) {
  const std::string exposure = e.get("exposure").as_string(),
                    outcome = e.get("outcome").as_string();
  const auto before = ancestors(g, g.index(outcome)),
             after = descendants(g, g.index(exposure));
  std::vector<int> result;
  const auto &edges = g.json.get("edges");
  for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
    const auto &edge = edges.at(i);
    const bool direct = edge.get("source").as_string() == exposure &&
                        edge.get("target").as_string() == outcome;
    if (direct)
      result.insert(result.begin(), i);
    else if (e.get("kind").as_string() != "coefficient" &&
             after[g.index(edge.get("source").as_string())] &&
             before[g.index(edge.get("target").as_string())])
      result.push_back(i);
  }
  return result;
}
J change_arrow(const Graph &g, int edge, double value) {
  J out = g.json, edges = J::array();
  for (int i = 0; i < static_cast<int>(g.json.get("edges").size()); ++i) {
    J item = g.json.get("edges").at(i);
    if (i == edge)
      item.set("coefficient", value);
    edges.push(item);
  }
  out.set("edges", edges);
  return out;
}
J null_graph(const Graph &g, const J &e) {
  const auto before = ancestors(g, g.index(e.get("outcome").as_string()));
  const std::string a = e.get("exposure").as_string();
  J out = g.json, edges = J::array(), nodes = J::array();
  for (J edge : g.json.get("edges").elements()) {
    if (edge.get("source").as_string() == a &&
        before[g.index(edge.get("target").as_string())])
      edge.set("coefficient", 0);
    edges.push(edge);
  }
  for (J node : g.json.get("nodes").elements()) {
    J terms = J::array();
    for (J term : node.get("interactions").elements()) {
      if (before[g.index(node.get("id").as_string())] &&
          (term.get("left").as_string() == a ||
           term.get("right").as_string() == a))
        term.set("coefficient", 0);
      terms.push(term);
    }
    node.set("interactions", terms);
    if (node.find("functions")) {
      J functions = J::array();
      for (J term : node.get("functions").elements()) {
        if (before[g.index(node.get("id").as_string())] &&
            (term.get("left").as_string() == a || term.get("right").as_string() == a))
          term.set("coefficient", 0);
        functions.push(term);
      }
      node.set("functions", functions);
    }
    if(node.find("lags")) {
      J lags=J::array(); for(J lag:node.get("lags").elements()) {
        if(before[g.index(node.get("id").as_string())] && lag.get("source").as_string()==a) lag.set("coefficient",0);
        lags.push(lag);
      } node.set("lags",lags);
    }
    nodes.push(node);
  }
  out.set("edges", edges);
  out.set("nodes", nodes);
  return out;
}
} // namespace ssc::dag

