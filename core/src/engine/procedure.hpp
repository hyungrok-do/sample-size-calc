// Procedure model: parameter schema, scenarios, result rows and the registry.
//
// A procedure declares its inputs (with declarative visibility rules that the
// UI evaluates too), the quantities it can solve for, and a `run` function that
// computes one fully specified scenario.  The engine (runner.cpp) takes care of
// parsing multi-valued inputs, expanding scenarios, running them in parallel,
// dropout inflation, plots and the JSON report.
#pragma once

#include <atomic>
#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

#include "util/parallel.hpp"

namespace ssc {

/// Error raised by a procedure for an infeasible scenario (reported in its row).
struct CalcError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// ------------------------------------------------------------------ visibility

struct Clause {
  std::string key;  // "$solve" or a select/bool parameter id
  std::vector<std::string> values;
  bool negate = false;
};
using Conj = std::vector<Clause>;  // all clauses must hold

/// Disjunction of conjunctions; empty = always visible.
struct Visibility {
  std::vector<Conj> any;
  bool always() const { return any.empty(); }
  template <class Lookup>
  bool eval(const Lookup& lookup) const {
    if (any.empty()) return true;
    for (const auto& conj : any) {
      bool ok = true;
      for (const auto& c : conj) {
        const std::string v = lookup(c.key);
        bool in = false;
        for (const auto& x : c.values) in = in || (x == v);
        if (in == c.negate) {
          ok = false;
          break;
        }
      }
      if (ok) return true;
    }
    return false;
  }
};

inline Clause solving(std::vector<std::string> v) { return Clause{"$solve", std::move(v), false}; }
inline Clause not_solving(std::vector<std::string> v) { return Clause{"$solve", std::move(v), true}; }
inline Clause is(std::string key, std::vector<std::string> v) { return Clause{std::move(key), std::move(v), false}; }
inline Clause is_not(std::string key, std::vector<std::string> v) { return Clause{std::move(key), std::move(v), true}; }
inline Visibility when(Conj c) { return Visibility{{std::move(c)}}; }
inline Visibility when(Clause c) { return Visibility{{Conj{std::move(c)}}}; }
inline Visibility when_any(std::vector<Conj> v) { return Visibility{std::move(v)}; }

// ------------------------------------------------------------------ parameters

enum class PType { Number, Select, Bool, Vector };

struct Option {
  std::string value;
  std::string label;
  Visibility visible;  // e.g. allocation choices that only apply when solving for N
};

struct ParamDef {
  std::string id;
  std::string label;
  std::string symbol;  // short header used in report tables (defaults to label)
  std::string group;
  std::string help;
  PType type = PType::Number;
  std::string def;  // default text ("true"/"false" for Bool)
  bool multi = true;
  bool integer = false;
  double min = -std::numeric_limits<double>::infinity();
  double max = std::numeric_limits<double>::infinity();
  bool min_open = false;
  bool max_open = false;
  std::vector<Option> options;
  std::vector<std::string> suggest;
  Visibility visible;

  // fluent setters
  ParamDef& grp(std::string g) { group = std::move(g); return *this; }
  ParamDef& sym(std::string s) { symbol = std::move(s); return *this; }
  ParamDef& hint(std::string h) { help = std::move(h); return *this; }
  ParamDef& range(double lo, double hi, bool lo_open = false, bool hi_open = false) {
    min = lo; max = hi; min_open = lo_open; max_open = hi_open; return *this;
  }
  ParamDef& above(double lo) { min = lo; min_open = true; return *this; }
  ParamDef& at_least(double lo) { min = lo; min_open = false; return *this; }
  ParamDef& whole() { integer = true; return *this; }
  ParamDef& single() { multi = false; return *this; }
  ParamDef& sugg(std::vector<std::string> s) { suggest = std::move(s); return *this; }
  ParamDef& show(Visibility v) { visible = std::move(v); return *this; }
};

ParamDef num(std::string id, std::string label, std::string def);
ParamDef select(std::string id, std::string label, std::vector<Option> options, std::string def);
ParamDef flag(std::string id, std::string label, bool def);
ParamDef vec(std::string id, std::string label, std::string def);

// ------------------------------------------------------------------ scenario

/// One fully specified set of inputs (one value per visible parameter).
class Scenario {
 public:
  std::string solve;
  std::unordered_map<std::string, double> num;
  std::unordered_map<std::string, std::string> sel;
  std::unordered_map<std::string, std::vector<double>> vec;

  double operator()(const std::string& id) const {
    auto it = num.find(id);
    if (it == num.end()) throw CalcError("missing input: " + id);
    return it->second;
  }
  double get(const std::string& id, double def) const {
    auto it = num.find(id);
    return it == num.end() ? def : it->second;
  }
  bool has(const std::string& id) const { return num.count(id) || sel.count(id) || vec.count(id); }
  const std::string& s(const std::string& id) const {
    auto it = sel.find(id);
    if (it == sel.end()) throw CalcError("missing selection: " + id);
    return it->second;
  }
  bool is(const std::string& id, const char* value) const {
    auto it = sel.find(id);
    return it != sel.end() && it->second == value;
  }
  bool flag(const std::string& id) const {
    auto it = sel.find(id);
    return it != sel.end() && it->second == "true";
  }
  const std::vector<double>& v(const std::string& id) const {
    auto it = vec.find(id);
    if (it == vec.end()) throw CalcError("missing list input: " + id);
    return it->second;
  }
};

// ------------------------------------------------------------------ results

using Cell = std::variant<std::monostate, double, std::string>;

struct Table {
  std::string title;
  std::vector<std::pair<std::string, std::string>> columns;  // key, label
  std::vector<std::string> formats;                          // per column
  std::vector<std::vector<Cell>> rows;
};

struct Row {
  std::vector<std::pair<std::string, Cell>> cells;
  std::string error;
  std::string detail;  // a numerical routine's own words behind a generic error (for tests; not shown)
  std::string statement;
  std::vector<std::string> notes;
  std::vector<Table> tables;

  void set(const std::string& key, double v);
  void set(const std::string& key, const std::string& v);
  void set(const std::string& key, const char* v) { set(key, std::string(v)); }
  double num(const std::string& key) const;  // NaN if absent / not numeric
  bool has(const std::string& key) const;
};

struct Column {
  std::string key;
  std::string label;
  std::string fmt;  // "prob" (power), "int", "num", "alpha", "pct", "text", "ratio"
  std::string tip;
};

struct SolveOption {
  std::string value;
  std::string label;
  std::string ykey;    // row key plotted as the result (defaults to value)
  std::string ylabel;  // axis label
};

/// Progress of one scenario in `units` steps of the token's total; simulation
/// procedures report it as they go and the runner completes it afterwards.
class Progress {
 public:
  Progress(Token* token, long units) : token_(token), units_(units) {}
  /// Report that `fraction` (0..1) of the scenario is done; never moves backwards.
  void set(double fraction) {
    if (!token_) return;
    const double f = fraction < 0 ? 0 : (fraction > 1 ? 1 : fraction);
    const long target = static_cast<long>(f * static_cast<double>(units_));
    long current = reported_.load(std::memory_order_relaxed);
    while (target > current)
      if (reported_.compare_exchange_weak(current, target, std::memory_order_relaxed)) {
        token_->add_done(target - current);
        return;
      }
  }

 private:
  Token* token_;
  long units_;
  std::atomic<long> reported_{0};
};

struct RunContext {
  const Token* token = nullptr;
  int threads = 1;
  bool approximate = false;     // simulation procedures: analytic approximation instead of simulating
  Progress* progress = nullptr;  // progress of the scenario being computed (may be null)
  int curve_points = 0;          // a point of the power curves, which take this many runs in all (0: not a curve)
  void check() const {
    if (token) token->check();
  }
  void report(double fraction) const {
    if (progress) progress->set(fraction);
  }
};

struct Procedure {
  std::string id;
  std::string name;
  std::string category;
  std::string subcategory;
  std::string summary;
  std::vector<std::string> keywords;
  std::vector<SolveOption> solve_for;
  std::vector<ParamDef> params;

  std::function<void(const Scenario&, Row&, RunContext&)> run;
  std::function<std::vector<Column>(const Scenario&)> columns;
  std::function<std::string(const Scenario&, const Row&)> statement;

  std::vector<std::string> references;
  std::vector<std::pair<std::string, std::string>> definitions;
  std::vector<std::string> notes;

  // power-curve support: vary `curve_param` (a sample-size input) with solve = "power"
  std::string curve_param;
  std::string curve_label;
  double curve_min = 2;
  bool curve_integer = true;

  // dropout inflation: group-size keys and total key in the result row
  std::vector<std::string> dropout_groups;
  std::string dropout_total;
  double dropout_multiple = 1;  // inflated sizes rounded up to a multiple of this (2: two sequences of equal size)

  std::string example_solve;
  std::map<std::string, std::string> example;
  int order = 100;  // ordering within a subcategory

  // simulation procedures estimate power from simulated studies; with
  // RunContext::approximate they return an analytic approximation instead
  bool simulation = false;

  // the time of one scenario relative to a typical one (unset: 1), when it can be much longer: the scenario
  // limit then applies to the sum, so that a request accepted is one that ends in reasonable time
  std::function<double(const Scenario&)> cost;
};

// ------------------------------------------------------------------ registry

std::vector<const Procedure*>& registry();
const Procedure* find_procedure(const std::string& id);

struct ProcRegistrar {
  explicit ProcRegistrar(Procedure (*make)());
};

/// The fixed order of navigator categories.
const std::vector<std::string>& category_order();

/// All procedures in catalog order (category, subcategory, order, name).
std::vector<const Procedure*> catalog_procedures();

}  // namespace ssc

#define SSC_CONCAT2(a, b) a##b
#define SSC_CONCAT(a, b) SSC_CONCAT2(a, b)
#define SSC_REGISTER(maker) static ::ssc::ProcRegistrar SSC_CONCAT(ssc_proc_registrar_, __LINE__)(maker)
