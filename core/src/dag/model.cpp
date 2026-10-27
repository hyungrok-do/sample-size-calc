#include "model.hpp"
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
      {"bernoulli", "binary"},     {"poisson", "count"},
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
    if (node.family == "gaussian" || node.family == "lognormal")
      node.sd = positive("sd");
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
Data simulate(const Graph &g, int n, std::uint64_t seed,
              const std::map<int, double> &interventions, bool observe,
              const Token *token) {
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
        eta += t.coefficient * data.values[t.left][i] *
               (t.right < 0 ? 1 : data.values[t.right][i]);
      if (!std::isfinite(eta))
        throw Error("generation_error",
                    node.label + ": a conditional predictor overflowed.");
      double value;
      if (node.family == "gaussian")
        value = eta + node.sd * rng.normal();
      else if (node.family == "uniform")
        value = eta + node.width * (rng.uniform() - .5);
      else if (node.family == "lognormal")
        value = std::exp(eta + node.sd * rng.normal());
      else if (node.family == "gamma")
        value = rng.gamma(node.shape) / node.shape * std::exp(eta);
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
      } else
        value = std::exp(std::log(node.scale) +
                         (std::log(rng.exponential()) - eta) / node.shape);
      if (!std::isfinite(value) ||
          ((node.family == "lognormal" || node.family == "gamma" ||
            node.kind == "survival") &&
           value <= 0))
        throw Error("generation_error",
                    node.label + ": a simulated value is outside the supported "
                                 "numerical range.");
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
  for (int j : g.order)
    order.push(g.nodes[j].id);
  return object({{"n", data.n},
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
  return d;
}
std::vector<bool> ancestors(const Graph &g, int node) {
  std::vector<bool> found(g.nodes.size());
  found[node] = true;
  for (auto it = g.order.rbegin(); it != g.order.rend(); ++it)
    if (found[*it])
      for (const auto &t : g.nodes[*it].terms) {
        found[t.left] = true;
        if (t.right >= 0)
          found[t.right] = true;
      }
  return found;
}
std::vector<bool> descendants(const Graph &g, int node) {
  std::vector<bool> found(g.nodes.size());
  found[node] = true;
  for (int j : g.order)
    for (const auto &t : g.nodes[j].terms)
      if (found[t.left] || (t.right >= 0 && found[t.right]))
        found[j] = true;
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
    nodes.push(node);
  }
  out.set("edges", edges);
  out.set("nodes", nodes);
  return out;
}
} // namespace ssc::dag
