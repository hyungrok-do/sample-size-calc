// The executable DAG uses the same random streams and thread pool as the
// calculator.
#pragma once
#include "util/json.hpp"
#include "util/parallel.hpp"
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace ssc::dag {
using J = json::Value;
using Vec = std::vector<double>;
using Columns = std::map<std::string, Vec>;
using Progress = std::function<void(const J &)>;
struct Error : std::runtime_error {
  std::string code;
  Error(std::string c, const std::string &message)
      : std::runtime_error(message), code(std::move(c)) {}
};
[[noreturn]] void fail(const std::string &message);
double number(const J &object, const std::string &key);
int integer(const J &object, const std::string &key, int lower, int upper,
            int fallback);
J object(std::initializer_list<std::pair<std::string, J>> values);
J array(const Vec &values);
struct Term {
  int left, right = -1;
  double coefficient = 0;
};
struct Node {
  std::string id, label, kind, family, censor;
  double intercept = 0, sd = 1, width = 1, shape = 1, scale = 1, cutoff = 20;
  int censor_index = -1;
  std::uint64_t stream = 0;
  std::vector<Term> terms;
};
struct Graph {
  J json;
  std::vector<Node> nodes;
  std::vector<int> order;
  int index(const std::string &id) const;
};
Graph parse_graph(const J &value);
struct Data {
  int n = 0;
  std::vector<Vec> values, times, events;
};
Data simulate(const Graph &graph, int n, std::uint64_t seed,
              const std::map<int, double> &interventions = {},
              bool observe = true, const Token *token = nullptr);
J data_json(const Graph &graph, const Data &data, std::uint64_t seed,
            bool summarize);
Data read_data(const Graph &graph, const J &value);
std::vector<bool> ancestors(const Graph &graph, int node);
std::vector<bool> descendants(const Graph &graph, int node);
std::vector<int> calibration_edges(const Graph &graph, const J &estimand);
J null_graph(const Graph &graph, const J &estimand);
J change_arrow(const Graph &graph, int edge, double value);
J run_request(const J &request, Token *token = nullptr,
              const Progress &progress = {});
} // namespace ssc::dag
