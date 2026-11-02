#pragma once
#include "statistics.hpp"
#include "sim/random.hpp"
namespace ssc::dag {
void parse_study(Graph &graph);
double draw_value(const Node &node, double eta, sim::Rng &rng, double normal = NAN);
Data simulate_study(const Graph &graph, int units, std::uint64_t seed,
                    const std::map<int,double> &interventions, bool observe, const Token *token);
void read_study_data(const Graph &graph, const J &value, Data &data);
J study_plan(const Graph &graph, const J &estimand, const J &analysis);
Vec study_weights(const Graph &graph, const J &estimand);
Fit analyze_study(const Graph &graph, const J &estimand, const Data &data, const J &analysis);
J study_population(const Graph &graph, const J &estimand, int units, std::uint64_t seed, const Token *token);
int study_rows(const Graph &graph, int units);
}
