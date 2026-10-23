#include "engine/runner.hpp"

#include "engine/search.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <new>
#include <set>

#include "engine/blocks.hpp"
#include "engine/parse.hpp"
#include "util/numtext.hpp"

namespace ssc {

namespace {

constexpr double MAX_SCENARIOS = 20000;  // combinations of input values in one request
constexpr std::size_t MAX_TABLE_ROWS = 1000;  // rows that keep their detail tables (the page draws 1,000 rows)

const char* type_name(PType t) {
  switch (t) {
    case PType::Number: return "number";
    case PType::Select: return "select";
    case PType::Bool: return "bool";
    case PType::Vector: return "vector";
  }
  return "number";
}

json::Value visibility_json(const Visibility& v) {
  json::Value any = json::Value::array();
  for (const auto& conj : v.any) {
    json::Value c = json::Value::array();
    for (const auto& cl : conj) {
      json::Value o = json::Value::object();
      o.set("key", cl.key);
      json::Value vals = json::Value::array();
      for (const auto& s : cl.values) vals.push(s);
      o.set("values", vals);
      o.set("negate", cl.negate);
      c.push(o);
    }
    any.push(c);
  }
  return any;
}

// The largest count (a group size, controls, events) an answer gives; past it no study could have it.
constexpr double MAX_COUNT = 1e15;

// The sizes an input that is not a count may have, other than 0: no study's quantity lies beyond them in any
// sensible unit, and within them the products, squares and ratios of the methods stay far from overflow and
// underflow.
constexpr double MIN_SIZE = 1e-30, MAX_SIZE = 1e30;

// The rule for a value outside those sizes, in the terms of the input's own limits: "0 or" only where its limits
// allow 0, "in size" where they allow negative values, and each end only where its limits reach past it.
std::string size_rule(const ParamDef& p) {
  const bool zero = (p.min_open ? p.min < 0 : p.min <= 0) && (p.max_open ? p.max > 0 : p.max >= 0);
  const bool small = p.min < MIN_SIZE && p.max > -MIN_SIZE, large = p.max > MAX_SIZE || p.min < -MAX_SIZE;
  std::string rule = small && large ? "between 1e-30 and 1e30" : (small ? "at least 1e-30" : "at most 1e30");
  if (p.min < 0) rule += " in size";
  return "values must be " + std::string(zero ? "0 or " : "") + rule;
}

json::Value number_or_null(double d) { return std::isfinite(d) ? json::Value(d) : json::Value(); }

json::Value cell_json(const Cell& c) {
  if (const double* d = std::get_if<double>(&c)) return number_or_null(*d);
  if (const std::string* s = std::get_if<std::string>(&c)) return json::Value(*s);
  return json::Value();
}

std::string normalize_bool(const std::string& s) {
  std::string t;
  for (char ch : s) t.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
  return (t == "true" || t == "1" || t == "yes" || t == "on") ? "true" : "false";
}

std::string value_label(double v) { return blocks::fnum(v, 6); }

std::string vector_label(const std::vector<double>& v) {
  std::string s;
  for (std::size_t i = 0; i < v.size(); ++i) {
    if (i) s += " ";
    s += value_label(v[i]);
  }
  return s;
}

struct Axis {
  const ParamDef* p = nullptr;
  std::vector<double> nums;
  std::vector<std::vector<double>> vecs;
  std::size_t size() const { return p->type == PType::Vector ? vecs.size() : nums.size(); }
};

std::string param_short(const ParamDef& p) { return p.symbol.empty() ? p.label : p.symbol; }

// Checks every test shares. With no effect at all a test rejects with probability alpha, so a target power at or
// below alpha needs no study; and a one-sided test at alpha of 0.5 or more rejects at least as often as not.
void check_levels(const Scenario& s) {
  const auto a = s.num.find("alpha");
  if (a == s.num.end()) return;
  const auto pw = s.num.find("power");
  if (s.solve != "power" && pw != s.num.end() && !(pw->second > a->second))
    throw CalcError("the target power must exceed alpha (with no effect at all, the test rejects with probability "
                    "alpha)");
  const auto side = s.sel.find("alternative");
  if (side != s.sel.end() && side->second != "two" && !(a->second < 0.5))
    throw CalcError("a one-sided test needs alpha below 0.5");
}

// Group sizes that follow from a ratio (N2 = R × N1, rounded up) or a percentage of a total are held to the
// smallest group a size typed in may have: a group of one is refused, as N2 = 1 typed in would be, and the
// search for the sample size (which starts where both groups are large enough) stays the smallest that works.
void check_groups(const Procedure& p, const Scenario& s) {
  const auto alloc = s.sel.find("alloc");
  if (alloc == s.sel.end()) return;
  const bool ratio = alloc->second == "ratio" && s.num.count("n1") && s.num.count("r");
  const bool pct = alloc->second == "pct" && s.num.count("ntotal") && s.num.count("pct1");
  if (!ratio && !pct) return;  // the sizes are typed in (checked against their limits) or solved for
  auto least = [&](const char* id) {
    for (const auto& q : p.params)
      if (q.id == id) return q.min;
    return 1.0;
  };
  const blocks::TwoN n = blocks::two_group_n(s);
  const double min2 = least("n2"), min1 = least("n1");
  if (n.n2 < min2)
    throw CalcError("group 2 would have N2 = " + blocks::fint(n.n2) +
                    (ratio ? " (R × N1, rounded up)" : " (the rest of the total)") + "; each group needs at least " +
                    blocks::fint(min2));
  if (n.n1 < min1)
    throw CalcError("group 1 would have N1 = " + blocks::fint(n.n1) + " (its percentage of the total); each group "
                    "needs at least " + blocks::fint(min1));
}

}  // namespace

// ---------------------------------------------------------------------------
json::Value procedure_json(const Procedure& p) {
  json::Value o = json::Value::object();
  o.set("id", p.id);
  o.set("name", p.name);
  o.set("category", p.category);
  o.set("subcategory", p.subcategory);
  o.set("summary", p.summary);
  json::Value kw = json::Value::array();
  for (const auto& k : p.keywords) kw.push(k);
  o.set("keywords", kw);
  json::Value sf = json::Value::array();
  for (const auto& s : p.solve_for) {
    json::Value x = json::Value::object();
    x.set("value", s.value);
    x.set("label", s.label);
    x.set("ykey", s.ykey.empty() ? s.value : s.ykey);  // the result row's key of the solved quantity
    sf.push(x);
  }
  o.set("solve_for", sf);
  o.set("simulation", p.simulation);
  json::Value groups = json::Value::array();
  std::set<std::string> seen;
  for (const auto& prm : p.params)
    if (seen.insert(prm.group).second) groups.push(prm.group);
  o.set("groups", groups);
  json::Value params = json::Value::array();
  for (const auto& prm : p.params) {
    json::Value x = json::Value::object();
    x.set("id", prm.id);
    x.set("label", prm.label);
    x.set("symbol", prm.symbol.empty() ? prm.label : prm.symbol);
    x.set("group", prm.group);
    x.set("help", prm.help);
    x.set("type", type_name(prm.type));
    x.set("default", prm.def);
    x.set("multi", prm.multi);
    x.set("integer", prm.integer);
    x.set("min", number_or_null(prm.min));
    x.set("max", number_or_null(prm.max));
    x.set("min_open", prm.min_open);
    x.set("max_open", prm.max_open);
    json::Value opts = json::Value::array();
    for (const auto& op : prm.options) {
      json::Value y = json::Value::object();
      y.set("value", op.value);
      y.set("label", op.label);
      y.set("visible", visibility_json(op.visible));
      opts.push(y);
    }
    x.set("options", opts);
    json::Value sug = json::Value::array();
    for (const auto& s : prm.suggest) sug.push(s);
    x.set("suggest", sug);
    x.set("visible", visibility_json(prm.visible));
    params.push(x);
  }
  o.set("params", params);
  json::Value refs = json::Value::array();
  for (const auto& r : p.references) refs.push(r);
  o.set("references", refs);
  json::Value defs = json::Value::array();
  for (const auto& d : p.definitions) {
    json::Value pair = json::Value::array();
    pair.push(d.first);
    pair.push(d.second);
    defs.push(pair);
  }
  o.set("definitions", defs);
  json::Value notes = json::Value::array();
  for (const auto& n : p.notes) notes.push(n);
  o.set("notes", notes);
  json::Value ex = json::Value::object();
  ex.set("solve_for", p.example_solve.empty() ? (p.solve_for.empty() ? "" : p.solve_for[0].value) : p.example_solve);
  json::Value exin = json::Value::object();
  for (const auto& kv : p.example) exin.set(kv.first, kv.second);
  ex.set("inputs", exin);
  o.set("example", ex);
  return o;
}

json::Value catalog_json() {
  const std::vector<const Procedure*> procs = catalog_procedures();
  const auto& cats = category_order();
  json::Value o = json::Value::object();
  json::Value catv = json::Value::array();
  std::set<std::string> present;
  for (const auto* p : procs) present.insert(p->category);
  for (const auto& c : cats)
    if (present.count(c)) catv.push(c);
  for (const auto& c : present)
    if (std::find(cats.begin(), cats.end(), c) == cats.end()) catv.push(c);
  o.set("categories", catv);
  json::Value list = json::Value::array();
  for (const auto* p : procs) {
    json::Value x = json::Value::object();
    x.set("id", p->id);
    x.set("name", p->name);
    x.set("category", p->category);
    x.set("subcategory", p->subcategory);
    x.set("summary", p->summary);
    json::Value kw = json::Value::array();
    for (const auto& k : p->keywords) kw.push(k);
    x.set("keywords", kw);
    x.set("simulation", p->simulation);
    list.push(x);
  }
  o.set("procedures", list);
  return o;
}

// ---------------------------------------------------------------------------
json::Value run_request(const json::Value& req, Token* token) {
  const auto t0 = std::chrono::steady_clock::now();
  json::Value out = json::Value::object();
  Token local_token;
  if (!token) token = &local_token;

  const std::string pid = req.get("procedure").as_string();
  const Procedure* proc = find_procedure(pid);
  if (!proc) {
    out.set("ok", false);
    out.set("error", "unknown procedure: " + parse::quoted(pid));
    return out;
  }
  out.set("procedure", pid);
  std::string solve = req.get("solve_for").to_text();
  if (solve.empty() && !proc->solve_for.empty()) solve = proc->solve_for[0].value;
  const SolveOption* solve_opt = nullptr;
  for (const auto& so : proc->solve_for)
    if (so.value == solve) solve_opt = &so;
  if (!solve_opt) {
    out.set("ok", false);
    out.set("error", "this procedure cannot solve for '" + parse::quoted(solve) + "'");
    return out;
  }
  out.set("solve_for", solve);

  const json::Value& inputs = req.get("inputs");
  const json::Value& options = req.get("options");
  // numbers from the request are clamped as numbers before they become integers (a cast of an out-of-range or
  // NaN double is undefined)
  auto option = [&](const char* key, double fallback, double lo, double hi) {
    const double v = options.get(key).as_number(fallback);
    return std::isfinite(v) ? std::clamp(std::floor(v), lo, hi) : fallback;
  };
  const int threads = resolve_threads(static_cast<int>(option("threads", 0, 0, 4096)));
  const bool want_curve = options.get("power_curve").as_bool(true);
  const int curve_points = static_cast<int>(option("curve_points", 41, 5, 400));
  // a caller can lower the limit of scenarios (the MCP server does), never raise it
  const std::size_t max_scen = static_cast<std::size_t>(option("max_scenarios", MAX_SCENARIOS, 1, MAX_SCENARIOS));
  const bool approximate = proc->simulation && options.get("approximate").as_bool(false);

  auto raw = [&](const ParamDef& p) -> std::string {
    const json::Value* v = inputs.find(p.id);
    if (v && !v->is_null()) {
      const std::string t = v->to_text();
      bool blank = true;
      for (char c : t)
        if (!std::isspace(static_cast<unsigned char>(c))) blank = false;
      if (!blank || p.type == PType::Select || p.type == PType::Bool) return t;
      return t;  // keep blank: reported as missing for visible numeric inputs
    }
    return p.def;
  };

  // --- choices (select / bool) used by visibility rules ----------------------
  std::unordered_map<std::string, std::string> choice;
  auto lookup = [&](const std::string& key) -> std::string {
    if (key == "$solve") return solve;
    auto it = choice.find(key);
    return it == choice.end() ? std::string() : it->second;
  };
  for (const auto& p : proc->params) {
    if (p.type == PType::Bool) choice[p.id] = normalize_bool(raw(p));
    if (p.type == PType::Select) choice[p.id] = raw(p);
  }
  // hidden select options fall back to the first visible option; a value that is no option at all is an
  // input error (rather than quietly the default)
  std::vector<std::pair<const ParamDef*, std::string>> unknown_choice;
  for (const auto& p : proc->params) {
    if (p.type != PType::Select) continue;
    const std::string cur = choice[p.id];
    const json::Value* given = inputs.find(p.id);
    bool known = false;
    for (const auto& op : p.options) known = known || op.value == cur;
    if (!known && given && !given->is_null() && cur.find_first_not_of(" \t\r\n") != std::string::npos)
      unknown_choice.emplace_back(&p, cur);
    bool ok = false;
    for (const auto& op : p.options)
      if (op.value == cur && op.visible.eval(lookup)) ok = true;
    if (!ok) {
      std::string fallback;
      for (const auto& op : p.options)
        if (op.value == p.def && op.visible.eval(lookup)) fallback = op.value;
      if (fallback.empty())
        for (const auto& op : p.options)
          if (op.visible.eval(lookup)) {
            fallback = op.value;
            break;
          }
      choice[p.id] = fallback;
    }
  }

  // --- parse visible inputs ---------------------------------------------------
  json::Value errors = json::Value::object();
  bool has_errors = false;
  for (const auto& u : unknown_choice) {
    if (!u.first->visible.eval(lookup)) continue;
    std::string list;
    for (const auto& op : u.first->options)
      if (op.visible.eval(lookup)) list += (list.empty() ? "" : ", ") + op.value;
    errors.set(u.first->id, "unknown choice '" + parse::quoted(u.second) + "': use one of " + list);
    has_errors = true;
  }
  // What the request gives that this calculation does not use, said rather than dropped (a program or a language
  // model that sends n2 with equal groups would otherwise believe it was used): names that are no input of the
  // procedure, inputs that the chosen options leave out, and choices not available here, replaced by the default.
  {
    json::Value unknown = json::Value::array(), ignored = json::Value::array(), replaced = json::Value::object();
    if (inputs.is_object())
      for (const auto& kv : inputs.items()) {
        const ParamDef* p = nullptr;
        for (const auto& q : proc->params)
          if (q.id == kv.first) p = &q;
        if (!p) {
          unknown.push(kv.first);
          continue;
        }
        const std::string given = kv.second.is_null() ? std::string() : kv.second.to_text();
        if (given.find_first_not_of(" \t\r\n") == std::string::npos) continue;
        if (!p->visible.eval(lookup)) {
          ignored.push(kv.first);
        } else if (p->type == PType::Select && given != choice[p->id] &&
                   std::any_of(p->options.begin(), p->options.end(),
                               [&](const Option& op) { return op.value == given; })) {
          json::Value r = json::Value::object();
          r.set("given", given);
          r.set("used", choice[p->id]);
          replaced.set(kv.first, r);
        }
      }
    if (unknown.size() > 0) out.set("unknown_inputs", unknown);
    if (ignored.size() > 0) out.set("ignored_inputs", ignored);
    if (replaced.size() > 0) out.set("replaced_choices", replaced);
  }
  Scenario base;
  base.solve = solve;
  std::vector<Axis> axes;
  for (const auto& p : proc->params) {
    if (!p.visible.eval(lookup)) continue;
    if (p.type == PType::Select || p.type == PType::Bool) {
      base.sel[p.id] = choice[p.id];
      continue;
    }
    const std::string text = raw(p);
    // other than 0, a number that is not a count lies within the sizes
    auto outside = [&](double v) {
      return !p.integer && v != 0 && !(std::fabs(v) >= MIN_SIZE && std::fabs(v) <= MAX_SIZE);
    };
    auto check_value = [&](double v) -> std::string {
      if (!std::isfinite(v)) return "values must be finite numbers";
      if (p.integer && v != std::floor(v)) return "values must be whole numbers";
      const auto bound = [](double b) { return blocks::fnum(b, 6); };
      if (p.min_open ? !(v > p.min) : !(v >= p.min))
        return std::string("values must be ") + (p.min_open ? "greater than " : "at least ") + bound(p.min);
      if (p.max_open ? !(v < p.max) : !(v <= p.max))
        return std::string("values must be ") + (p.max_open ? "less than " : "at most ") + bound(p.max);
      if (outside(v)) return size_rule(p);
      // a count (subjects, clusters, levels ...) up to a billion stays exact in the integer arithmetic of the
      // exact methods; no study needs more
      if (p.integer && v > 1e9) return "values must be at most 1,000,000,000";
      return {};
    };
    if (p.type == PType::Number) {
      parse::Numbers r = parse::numbers(text);
      std::string err = r.error;
      if (err.empty() && !p.multi && r.values.size() > 1) err = "only one value is allowed";
      if (err.empty())
        for (double v : r.values) {
          err = check_value(v);
          if (!err.empty()) break;
        }
      // so does each end of a range as typed, even one the steps do not reach ('0 to 1.5e30 by 1e30' stops at
      // 1e30), though only that rule: an input's own limits hold for the values (power '0.5 to 1 by 0.2')
      for (double v : r.ends)
        if (err.empty() && outside(v)) err = size_rule(p);
      if (!err.empty()) {
        errors.set(p.id, err);
        has_errors = true;
        continue;
      }
      if (r.values.size() == 1)
        base.num[p.id] = r.values[0];
      else
        axes.push_back(Axis{&p, r.values, {}});
    } else {  // vector
      parse::Vectors r = parse::vectors(text);
      std::string err = r.error;
      if (err.empty())
        for (const auto& vv : r.vectors)
          for (double v : vv)
            if (err.empty()) err = check_value(v);
      for (double v : r.ends)
        if (err.empty() && outside(v)) err = size_rule(p);
      if (!err.empty()) {
        errors.set(p.id, err);
        has_errors = true;
        continue;
      }
      if (r.vectors.size() == 1)
        base.vec[p.id] = r.vectors[0];
      else
        axes.push_back(Axis{&p, {}, r.vectors});
    }
  }
  if (has_errors) {
    out.set("ok", false);
    out.set("error", "Some inputs are invalid.");
    out.set("input_errors", errors);
    return out;
  }

  // --- scenario expansion (last axis varies fastest) -------------------------
  std::size_t count = 1;
  for (const auto& a : axes) {
    count *= a.size();
    if (count > max_scen) {
      out.set("ok", false);
      out.set("error", "Too many combinations of input values (limit " + std::to_string(max_scen) + ").");
      return out;
    }
  }
  // each scenario holds its own copy of the lists, and a row shows them: bound their values over all scenarios
  std::size_t listed = 0;
  for (const auto& kv : base.vec) listed += kv.second.size();
  for (const auto& a : axes) {
    std::size_t longest = 0;
    for (const auto& v : a.vecs) longest = std::max(longest, v.size());
    listed += longest;
  }
  if (count > 1 && listed * count > 1000000) {
    out.set("ok", false);
    out.set("error", "Too many values in the lists over all the combinations of input values (limit 1,000,000): use "
                     "fewer combinations or shorter lists.");
    return out;
  }
  std::vector<Scenario> scen(count, base);
  std::vector<std::vector<std::size_t>> index(count, std::vector<std::size_t>(axes.size()));
  for (std::size_t i = 0; i < count; ++i) {
    std::size_t rem = i;
    for (std::size_t k = axes.size(); k-- > 0;) {
      const std::size_t sz = axes[k].size();
      const std::size_t j = rem % sz;
      rem /= sz;
      index[i][k] = j;
      if (axes[k].p->type == PType::Vector)
        scen[i].vec[axes[k].p->id] = axes[k].vecs[j];
      else
        scen[i].num[axes[k].p->id] = axes[k].nums[j];
    }
  }

  if (proc->cost && count > 1) {  // one scenario is always computed (the count limit has let it through)
    const double cost = proc->cost(scen);
    if (cost > static_cast<double>(max_scen)) {
      out.set("ok", false);
      out.set("error", "Too many combinations of input values for designs this long to compute (limit about " +
                           std::to_string(std::max(1L, static_cast<long>(static_cast<double>(max_scen) * count / cost))) +
                           " here).");
      return out;
    }
  }

  // --- compute rows in parallel -------------------------------------------
  const bool curve_ok = want_curve && !proc->curve_param.empty();
  const std::size_t n_curves = curve_ok ? std::min<std::size_t>(count, 8) : 0;
  const long units = proc->simulation && !approximate ? 1000 : 1;  // progress steps per scenario
  token->set_total(static_cast<long>(count) * units + static_cast<long>(n_curves) * curve_points);
  std::vector<Row> rows(count);
  std::atomic<bool> tables_dropped{false}, curve_approximated{false}, curve_out_of_memory{false};
  bool curves_drawn = false;
  RunContext ctx{token, threads, approximate, nullptr};
  auto answer_key = [&](const std::string& solve) {
    for (const auto& o : proc->solve_for)
      if (o.value == solve) return o.ykey.empty() ? o.value : o.ykey;
    return solve;
  };
  auto cell = [](const Row& r, const std::string& key) -> const double* {
    for (const auto& kv : r.cells)
      if (kv.first == key) return std::get_if<double>(&kv.second);
    return nullptr;
  };
  auto nan_cell = [&](const Row& r, const std::string& key) {  // no number, or an infinite one
    const double* v = cell(r, key);
    return v && !std::isfinite(*v);
  };
  // the counts of the table (group sizes, controls, clusters): one that overflowed, or that no study could have
  // (past MAX_COUNT; past 2^53 a double does not even hold whole numbers), says so rather than being shown
  std::vector<std::string> count_keys;
  if (proc->columns && count > 0)
    for (const auto& c : proc->columns(scen[0]))
      if (c.fmt == "int") count_keys.push_back(c.key);
  auto bad_count = [&](const Row& r, const std::string& key) {
    const double* v = cell(r, key);
    if (!v) return false;
    if (std::isnan(*v)) return key == "n" || key == "n1" || key == "n2";  // a criterion that does not apply is blank
    return std::fabs(*v) > MAX_COUNT;
  };
  try {
    parallel_for(
        count, threads,
        [&](std::size_t i) {
          Row r;
          Progress progress(token, units);
          RunContext local = ctx;
          local.progress = &progress;
          try {
            check_levels(scen[i]);
            check_groups(*proc, scen[i]);
            proc->run(scen[i], r, local);
            // a probability computed as 1 + a rounding error is shown as 1
            const double pw = r.num("power");
            if (pw > 1 && pw < 1 + 1e-9) r.set("power", 1.0);
            if (pw < 0 && pw > -1e-9) r.set("power", 0.0);
            // an answer that came out as no number or an infinite one (an overflow or an undefined quantity at
            // extreme inputs) is said, rather than shown as a blank
            for (const std::string& key : {answer_key(scen[i].solve), std::string("power")})
              if (nan_cell(r, key))
                throw CalcError("the answer cannot be computed for these inputs; check for an extreme value");
            for (const std::string& key : count_keys)
              if (bad_count(r, key))
                throw CalcError("the answer cannot be computed for these inputs: a count would be infinite or "
                                "beyond 10^15; check for an extreme value");
            if (proc->statement) r.statement = proc->statement(scen[i], r);
            // 20,000 rows of group-sequential boundaries would take a gigabyte
            if (i >= MAX_TABLE_ROWS && !r.tables.empty()) {
              r.tables.clear();
              tables_dropped = true;
            }
          } catch (const CalcError& e) {
            r.error = e.what();
          } catch (const Cancelled&) {
            throw;
          } catch (const std::bad_alloc&) {
            // a study, table or enumeration larger than the memory left (the page's engine has at most 2 GB): said
            // as it is, not as an extreme value
            r.error = "there is not enough memory to compute the answer for these inputs";
          } catch (const std::exception& e) {
            // a numerical routine gave up (a root it could not bracket, a series that did not converge): at the
            // inputs' extremes, as the stress sweeps show
            r.error = "the answer cannot be computed for these inputs; check for an extreme value";
            r.detail = e.what();
          }
          rows[i] = std::move(r);
          progress.set(1);
        },
        token);
  } catch (const Cancelled&) {
    out.set("ok", false);
    out.set("cancelled", true);
    out.set("error", "Calculation cancelled.");
    return out;
  }

  // building a large answer takes time too: Stop ends it
  auto stopped = [&] {
    out.set("ok", false);
    out.set("cancelled", true);
    out.set("error", "Calculation cancelled.");
    return out;
  };

  // --- columns ---------------------------------------------------------------
  std::vector<Column> columns = proc->columns ? proc->columns(scen[0]) : std::vector<Column>{};
  // columns without a value in any row (such as the Monte Carlo columns of an approximation) are left out
  if (std::any_of(rows.begin(), rows.end(), [](const Row& r) { return r.error.empty(); }))
    columns.erase(std::remove_if(columns.begin(), columns.end(),
                                 [&](const Column& c) {
                                   return std::none_of(rows.begin(), rows.end(),
                                                       [&](const Row& r) { return r.has(c.key); });
                                 }),
                  columns.end());
  // an input given several values that no column shows (the loss to follow-up, the cluster ratio ...): a column
  // of its own, so that the scenarios can be told apart
  for (std::size_t k = 0; k < axes.size(); ++k) {
    const ParamDef& p = *axes[k].p;
    if (std::any_of(columns.begin(), columns.end(),
                    [&](const Column& c) { return c.key == p.id || c.key == "target_" + p.id; }))
      continue;
    const bool vec = p.type == PType::Vector;
    columns.push_back({p.id, param_short(p), vec ? "text" : (p.id == "dropout" ? "pct" : (p.integer ? "int" : "num")),
                       p.label});
    for (std::size_t i = 0; i < count; ++i) {
      if (i % 64 == 0 && token->cancelled()) return stopped();
      if (rows[i].has(p.id)) continue;
      if (vec)
        rows[i].set(p.id, vector_label(axes[k].vecs[index[i][k]]));
      else
        rows[i].set(p.id, axes[k].nums[index[i][k]]);
    }
  }
  json::Value cols = json::Value::array();
  for (const auto& c : columns) {
    json::Value x = json::Value::object();
    x.set("key", c.key);
    x.set("label", c.label);
    x.set("fmt", c.fmt);
    x.set("tip", c.tip);
    cols.push(x);
  }
  out.set("columns", cols);

  // --- rows ------------------------------------------------------------------
  json::Value rowsj = json::Value::array();
  for (std::size_t i = 0; i < count; ++i) {
    if (i % 64 == 0 && token->cancelled()) return stopped();
    const Row& r = rows[i];
    json::Value x = json::Value::object();
    for (const auto& kv : r.cells) x.set(kv.first, cell_json(kv.second));
    if (!r.error.empty()) x.set("_error", r.error);
    if (!r.detail.empty()) x.set("_detail", r.detail);
    if (!r.statement.empty()) x.set("_statement", r.statement);
    if (!r.notes.empty()) {
      json::Value nn = json::Value::array();
      for (const auto& n : r.notes) nn.push(n);
      x.set("_notes", std::move(nn));
    }
    if (!r.tables.empty()) {
      json::Value tabs = json::Value::array();
      for (const auto& t : r.tables) {
        json::Value tj = json::Value::object();
        tj.set("title", t.title);
        json::Value tc = json::Value::array();
        for (std::size_t k = 0; k < t.columns.size(); ++k) {
          json::Value c = json::Value::object();
          c.set("key", t.columns[k].first);
          c.set("label", t.columns[k].second);
          c.set("fmt", k < t.formats.size() ? t.formats[k] : std::string("num"));
          tc.push(c);
        }
        tj.set("columns", std::move(tc));
        json::Value tr = json::Value::array();
        for (const auto& row : t.rows) {
          json::Value rr = json::Value::array();
          for (const auto& cell : row) rr.push(cell_json(cell));
          tr.push(std::move(rr));
        }
        tj.set("rows", std::move(tr));
        tabs.push(std::move(tj));
      }
      x.set("_tables", std::move(tabs));
    }
    rowsj.push(std::move(x));
  }
  out.set("rows", std::move(rowsj));

  // --- dropout-inflated sample sizes ------------------------------------------
  if (!proc->dropout_groups.empty()) {
    auto label_of = [&](const std::string& key) {
      for (const auto& c : columns)
        if (c.key == key) return c.label;
      return key;
    };
    json::Value dt = json::Value::object();
    json::Value dcols = json::Value::array();
    auto addc = [&](const std::string& key, const std::string& label, const std::string& fmt) {
      json::Value c = json::Value::object();
      c.set("key", key);
      c.set("label", label);
      c.set("fmt", fmt);
      dcols.push(c);
    };
    addc("dropout", "Dropout Rate", "pct");
    std::vector<std::string> keys = proc->dropout_groups;
    const bool has_total = !proc->dropout_total.empty() && proc->dropout_groups.size() > 1;
    for (const auto& k : keys) addc(k, label_of(k), "int");
    if (has_total) addc(proc->dropout_total, label_of(proc->dropout_total), "int");
    for (const auto& k : keys) addc(k + "_inflated", label_of(k) + "'", "int");
    if (has_total) addc(proc->dropout_total + "_inflated", label_of(proc->dropout_total) + "'", "int");
    auto dlabel = [&](const std::string& k) {
      const std::string l = label_of(k);
      return (!l.empty() && l[0] == 'N') ? "D" + l.substr(1) : "D(" + l + ")";
    };
    for (const auto& k : keys) addc(k + "_dropouts", dlabel(k), "int");
    if (has_total) addc(proc->dropout_total + "_dropouts", dlabel(proc->dropout_total), "int");
    dt.set("columns", dcols);
    json::Value drows = json::Value::array();
    for (std::size_t i = 0; i < count; ++i) {
      const double d = scen[i].get("dropout", 0);
      if (!(d > 0) || !rows[i].error.empty()) continue;
      json::Value x = json::Value::object();
      x.set("_row", static_cast<double>(i));
      x.set("dropout", d);
      double tot = 0, tot_inf = 0;
      bool ok = true;
      for (const auto& k : keys) {
        const double n = rows[i].num(k);
        if (!std::isfinite(n)) ok = false;
        const double mult = proc->dropout_multiple;
        const double ni = std::ceil(n / (1 - d) / mult - 1e-9) * mult;
        if (!(ni <= MAX_COUNT)) ok = false;  // a dropout rate near 1 inflates past any study: no row, as for the table
        x.set(k, n);
        x.set(k + "_inflated", ni);
        x.set(k + "_dropouts", ni - n);
        tot += n;
        tot_inf += ni;
      }
      if (!ok || (has_total && !(tot_inf <= MAX_COUNT))) continue;
      if (has_total) {
        x.set(proc->dropout_total, tot);
        x.set(proc->dropout_total + "_inflated", tot_inf);
        x.set(proc->dropout_total + "_dropouts", tot_inf - tot);
      }
      drows.push(x);
    }
    dt.set("rows", drows);
    if (drows.size() > 0) out.set("dropout", dt);
  }

  // --- plots -----------------------------------------------------------------
  json::Value plots = json::Value::array();
  const std::string ykey = solve_opt->ykey.empty() ? solve_opt->value : solve_opt->ykey;
  const std::string ylabel = solve_opt->ylabel.empty() ? solve_opt->label : solve_opt->ylabel;
  auto axis_value = [&](std::size_t k, std::size_t i) {
    return param_short(*axes[k].p) + " = " +
           (axes[k].p->type == PType::Vector ? vector_label(axes[k].vecs[index[i][k]])
                                              : value_label(axes[k].nums[index[i][k]]));
  };

  // Result plots: x = the numeric input with the most values (ties: the later
  // one, usually the effect size), series = the remaining input with the most
  // values; any further varying inputs produce separate panels.
  std::vector<std::size_t> numeric_axes;
  for (std::size_t k = 0; k < axes.size(); ++k)
    if (axes[k].p->type == PType::Number) numeric_axes.push_back(k);
  if (!numeric_axes.empty()) {
    std::size_t xa = numeric_axes[0];
    for (std::size_t k : numeric_axes)
      if (axes[k].size() >= axes[xa].size()) xa = k;
    std::size_t sa = axes.size();
    for (std::size_t k = 0; k < axes.size(); ++k)
      if (k != xa && (sa == axes.size() || axes[k].size() > axes[sa].size())) sa = k;
    std::map<std::vector<std::size_t>, std::vector<std::size_t>> panels;
    for (std::size_t i = 0; i < count; ++i) {
      std::vector<std::size_t> key;
      for (std::size_t k = 0; k < axes.size(); ++k)
        if (k != xa && k != sa) key.push_back(index[i][k]);
      panels[key].push_back(i);
    }
    std::size_t panel_no = 0;
    for (const auto& kv : panels) {
      if (++panel_no > 12) break;
      json::Value pl = json::Value::object();
      pl.set("kind", "result");
      pl.set("title", ylabel + " vs " + param_short(*axes[xa].p));
      std::string sub;
      const std::size_t i0 = kv.second[0];
      for (std::size_t k = 0; k < axes.size(); ++k) {
        if (k == xa || k == sa) continue;
        if (!sub.empty()) sub += ", ";
        sub += axis_value(k, i0);
      }
      pl.set("subtitle", sub);
      pl.set("x_label", axes[xa].p->label);
      pl.set("y_label", ylabel);
      pl.set("x_key", axes[xa].p->id);
      pl.set("y_key", ykey);
      std::map<std::size_t, std::vector<std::size_t>> series;
      for (std::size_t i : kv.second) series[sa < axes.size() ? index[i][sa] : 0].push_back(i);
      json::Value sj = json::Value::array();
      for (const auto& srs : series) {
        json::Value one = json::Value::object();
        one.set("name", sa < axes.size() ? axis_value(sa, srs.second[0]) : std::string());
        json::Value xs = json::Value::array(), ys = json::Value::array();
        for (std::size_t i : srs.second) {
          const double yv = rows[i].num(ykey);
          if (!rows[i].error.empty() || !std::isfinite(yv)) continue;
          xs.push(axes[xa].nums[index[i][xa]]);
          ys.push(yv);
        }
        one.set("x", xs);
        one.set("y", ys);
        sj.push(one);
      }
      pl.set("series", sj);
      plots.push(pl);
    }
  }

  // Power curves: power as a function of the sample size for (up to 8)
  // distinct designs.  Scenarios that differ only in the target power share
  // one curve, which then carries one marker per solution.
  if (curve_ok) {
    struct CurveJob {
      Scenario sc, last;  // its scenario, and that of its largest point
      std::string param, key, name;
      std::vector<double> xs, ys;
      std::vector<std::pair<double, double>> markers;
    };
    std::vector<CurveJob> jobs;
    std::string curve_label = proc->curve_label;
    std::string x_with_alloc;  // the allocation a click on the curve must also choose (N1 or N2 given alone)
    const bool solving_effect = solve != "n" && solve != "power";
    for (std::size_t i = 0; i < count && jobs.size() < n_curves; ++i) {
      if (!rows[i].error.empty()) continue;
      CurveJob job;
      job.sc = scen[i];
      job.sc.solve = "power";
      job.param = proc->curve_param;
      // a choice that the solved effect hides (its scale, say) takes its default in the power runs
      auto power_lookup = [&](const std::string& k) -> std::string {
        if (k == "$solve") return "power";
        auto it = job.sc.sel.find(k);
        return it == job.sc.sel.end() ? std::string() : it->second;
      };
      for (const auto& p : proc->params)
        if ((p.type == PType::Select || p.type == PType::Bool) && !job.sc.sel.count(p.id) &&
            p.visible.eval(power_lookup))
          job.sc.sel[p.id] = p.type == PType::Bool ? normalize_bool(p.def) : p.def;
      for (const auto& kv : rows[i].cells) {  // carry solved quantities back into the inputs
        const double* d = std::get_if<double>(&kv.second);
        if (!d || !std::isfinite(*d)) continue;
        for (const auto& p : proc->params)
          if (p.id == kv.first && p.type == PType::Number && !job.sc.num.count(p.id)) job.sc.num[p.id] = *d;
      }
      if (solving_effect && proc->curve_inputs) proc->curve_inputs(job.sc, rows[i]);
      auto alloc = job.sc.sel.find("alloc");
      if (alloc != job.sc.sel.end()) {
        const std::string a = alloc->second;
        if (a == "n1") {  // N1 fixed and N2 solved: vary N2
          alloc->second = "sep";
          job.sc.num["n1"] = rows[i].num("n1");
          job.sc.num["n2"] = rows[i].num("n2");
          job.param = "n2";
          curve_label = "N2 (Group 2 Sample Size)";
          x_with_alloc = "sep";
        } else if (a == "n2") {
          alloc->second = "sep";
          job.sc.num["n1"] = rows[i].num("n1");
          job.sc.num["n2"] = rows[i].num("n2");
          x_with_alloc = "sep";
        } else if (a == "pct") {
          job.sc.num["ntotal"] = rows[i].num("n");
          job.param = "ntotal";
          curve_label = "N (Total Sample Size)";
        }
      }
      const double c = job.sc.get(job.param, NAN);
      if (!std::isfinite(c)) continue;
      // identity of the curve: the original inputs except the target power and
      // the plotted sample size (+ the solved effect when solving an effect)
      std::string key;
      std::map<std::string, double> nums(scen[i].num.begin(), scen[i].num.end());
      for (const auto& kv : nums)
        if (kv.first != job.param && !(solve == "n" && kv.first == "power"))
          key += kv.first + "=" + text::format_roundtrip(kv.second) + ";";
      std::map<std::string, std::string> sels(scen[i].sel.begin(), scen[i].sel.end());
      for (const auto& kv : sels) key += kv.first + "=" + kv.second + ";";
      std::map<std::string, std::vector<double>> vecs(scen[i].vec.begin(), scen[i].vec.end());
      for (const auto& kv : vecs) key += kv.first + "=" + vector_label(kv.second) + ";";
      if (solving_effect) key += "effect=" + text::format_roundtrip(rows[i].num(ykey)) + ";";
      const double pw = rows[i].num("power");
      bool merged = false;
      for (auto& other : jobs)
        if (other.key == key) {
          other.markers.emplace_back(c, pw);
          merged = true;
        }
      if (merged) continue;
      std::string name;
      for (std::size_t k = 0; k < axes.size(); ++k) {
        const std::string& id = axes[k].p->id;
        if (id == job.param || (solve == "n" && id == "power")) continue;
        if (!name.empty()) name += ", ";
        name += axis_value(k, i);
      }
      if (solving_effect) name += std::string(name.empty() ? "" : ", ") + ylabel + " = " + value_label(rows[i].num(ykey));
      job.name = name;
      job.key = key;
      job.markers.emplace_back(c, pw);
      jobs.push_back(std::move(job));
    }
    for (auto& job : jobs) {
      double cmax = 0;
      for (const auto& m : job.markers) cmax = std::max(cmax, m.first);
      double lo = proc->curve_min;  // and no smaller than the input allows
      for (const auto& p : proc->params)
        if (p.id == job.param && std::isfinite(p.min))
          lo = std::max(lo, p.min_open ? (p.integer ? std::floor(p.min) + 1 : std::nextafter(p.min, INFINITY)) : p.min);
      double hi = std::max(2 * cmax, lo + 10);
      for (const auto& p : proc->params)  // and no larger than the input allows (whole numbers: up to 1e9)
        if (p.id == job.param) hi = std::max(lo, std::min({hi, p.max, p.integer ? 1e9 : INFINITY}));
      std::set<double> grid;
      for (int k = 0; k < curve_points; ++k) {
        double x = lo + (hi - lo) * k / (curve_points - 1);
        if (proc->curve_integer) x = std::round(x);
        grid.insert(x);
      }
      for (const auto& m : job.markers) grid.insert(m.first);
      job.xs.assign(grid.begin(), grid.end());
      job.ys.assign(job.xs.size(), NAN);
      job.last = job.sc;
      job.last.num[job.param] = job.xs.back();
    }
    std::vector<std::pair<std::size_t, std::size_t>> work;
    for (std::size_t j = 0; j < jobs.size(); ++j)
      for (std::size_t k = 0; k < jobs[j].xs.size(); ++k) work.emplace_back(j, k);
    // the curves of simulation procedures come from their analytic approximation
    RunContext curve_ctx = ctx;
    curve_ctx.approximate = proc->simulation;
    curve_ctx.curve_points = static_cast<int>(work.size());
    curve_ctx.curve_approximated = &curve_approximated;
    try {
      parallel_for(
          work.size(), threads,
          [&](std::size_t w) {
            CurveJob& job = jobs[work[w].first];
            Scenario sc = job.sc;
            sc.num[job.param] = job.xs[work[w].second];
            RunContext point = curve_ctx;  // a series takes one method throughout: the one its largest point needs
            point.curve_last = &job.last;
            Row r;
            try {
              proc->run(sc, r, point);
              job.ys[work[w].second] = r.num("power");
            } catch (const CalcError&) {
            } catch (const Cancelled&) {
              throw;
            } catch (const std::bad_alloc&) {
              curve_out_of_memory = true;  // a point left out, as any that cannot be computed, but said
            } catch (const std::exception&) {
            }
            token->add_done();
          },
          token);
    } catch (const Cancelled&) {
      out.set("ok", false);
      out.set("cancelled", true);
      out.set("error", "Calculation cancelled.");
      return out;
    }
    // each answer is marked at its own power: on its curve when the curve is computed as the answer is, near it
    // when one of the two is an approximation (the notes say which); a curve takes one method throughout
    if (!jobs.empty()) {
      curves_drawn = true;
      json::Value pl = json::Value::object();
      pl.set("kind", "curve");
      pl.set("x_key", jobs[0].param);  // the input along the x axis
      json::Value with = json::Value::object();  // and the choice to make with it
      if (!x_with_alloc.empty()) with.set("alloc", x_with_alloc);
      pl.set("x_with", with);
      pl.set("title", "Power vs " + curve_label);
      pl.set("subtitle", jobs.size() == 1 ? jobs[0].name : std::string());
      pl.set("x_label", curve_label);
      pl.set("y_label", "Power");
      json::Value sj = json::Value::array();
      for (const auto& job : jobs) {
        json::Value one = json::Value::object();
        one.set("name", job.name.empty() ? std::string("Power") : job.name);
        json::Value xs = json::Value::array(), ys = json::Value::array();
        for (std::size_t k = 0; k < job.xs.size(); ++k)
          if (std::isfinite(job.ys[k])) {
            xs.push(job.xs[k]);
            ys.push(job.ys[k]);
          }
        one.set("x", xs);
        one.set("y", ys);
        json::Value mks = json::Value::array();
        for (const auto& m : job.markers)
          if (std::isfinite(m.second)) {
            json::Value mk = json::Value::array();
            mk.push(m.first);
            mk.push(m.second);
            mks.push(mk);
          }
        one.set("markers", mks);
        sj.push(one);
      }
      pl.set("series", sj);
      plots.push(pl);
    }
  }
  out.set("plots", plots);

  // --- static report material ---------------------------------------------------
  json::Value refs = json::Value::array();
  for (const auto& r : proc->references) refs.push(r);
  out.set("references", refs);
  json::Value defs = json::Value::array();
  for (const auto& d : proc->definitions) {
    json::Value pr = json::Value::array();
    pr.push(d.first);
    pr.push(d.second);
    defs.push(pr);
  }
  out.set("definitions", defs);
  json::Value notes = json::Value::array();
  for (const auto& n : proc->notes) notes.push(n);
  if (curve_approximated)
    notes.push("Where enumerating every point of a power curve would take long, the curve uses the normal "
               "approximation; the answers marked on it are those of the table.");
  if (proc->simulation && !approximate && curves_drawn)
    notes.push("The power curve comes from the analytic approximation; the simulated answers are marked on it.");
  if (curve_out_of_memory)
    notes.push("Some points of the power curve are left out: there was not enough memory to compute them.");
  if (tables_dropped) {
    notes.push("Detail tables are given for the first 1,000 scenarios.");  // MAX_TABLE_ROWS
    out.set("tables_rows", static_cast<double>(MAX_TABLE_ROWS));
  }
  out.set("notes", notes);
  out.set("name", proc->name);
  out.set("solve_label", solve_opt->label);
  out.set("scenario_count", static_cast<double>(count));
  out.set("simulation", proc->simulation);
  out.set("approximate", approximate);
  out.set("threads", threads);
  out.set("elapsed_ms",
          std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
  out.set("ok", true);
  return out;
}

json::Value search_json(const std::string& query, std::size_t limit) {
  const search::Result res = search::query(query, limit);
  json::Value o = json::Value::object();
  o.set("query", res.query);
  o.set("corrected", res.corrected);
  o.set("partial", res.partial);
  json::Value hits = json::Value::array();
  for (const auto& h : res.hits) {
    const Procedure* p = find_procedure(h.id);
    json::Value x = json::Value::object();
    x.set("id", h.id);
    x.set("name", p ? p->name : h.id);
    x.set("category", p ? p->category : "");
    x.set("score", h.score);
    json::Value terms = json::Value::array();
    for (const auto& t : h.terms) terms.push(t);
    x.set("terms", terms);
    hits.push(x);
  }
  o.set("hits", hits);
  return o;
}

}  // namespace ssc
