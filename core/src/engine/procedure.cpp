#include "engine/procedure.hpp"

#include <algorithm>

#include <cmath>

namespace ssc {

ParamDef num(std::string id, std::string label, std::string def) {
  ParamDef p;
  p.id = std::move(id);
  p.label = std::move(label);
  p.def = std::move(def);
  p.type = PType::Number;
  return p;
}

ParamDef select(std::string id, std::string label, std::vector<Option> options, std::string def) {
  ParamDef p;
  p.id = std::move(id);
  p.label = std::move(label);
  p.type = PType::Select;
  p.options = std::move(options);
  p.def = std::move(def);
  p.multi = false;
  return p;
}

ParamDef flag(std::string id, std::string label, bool def) {
  ParamDef p;
  p.id = std::move(id);
  p.label = std::move(label);
  p.type = PType::Bool;
  p.def = def ? "true" : "false";
  p.multi = false;
  return p;
}

ParamDef vec(std::string id, std::string label, std::string def) {
  ParamDef p;
  p.id = std::move(id);
  p.label = std::move(label);
  p.type = PType::Vector;
  p.def = std::move(def);
  return p;
}

void Row::set(const std::string& key, double v) {
  for (auto& kv : cells)
    if (kv.first == key) {
      kv.second = v;
      return;
    }
  cells.emplace_back(key, v);
}

void Row::set(const std::string& key, const std::string& v) {
  for (auto& kv : cells)
    if (kv.first == key) {
      kv.second = v;
      return;
    }
  cells.emplace_back(key, v);
}

double Row::num(const std::string& key) const {
  for (const auto& kv : cells)
    if (kv.first == key) {
      if (const double* d = std::get_if<double>(&kv.second)) return *d;
      return std::nan("");
    }
  return std::nan("");
}

bool Row::has(const std::string& key) const {
  for (const auto& kv : cells)
    if (kv.first == key) return true;
  return false;
}

std::vector<const Procedure*>& registry() {
  // never destroyed, like the procedures it holds (see ProcRegistrar), so they stay reachable to the end and a
  // leak checker reports only real leaks
  static auto* procs = new std::vector<const Procedure*>;
  return *procs;
}

ProcRegistrar::ProcRegistrar(Procedure (*make)()) {
  // Procedures live for the whole program; leaking them avoids static
  // destruction-order issues when the library is unloaded.
  auto* p = new Procedure(make());
  registry().push_back(p);
}

const Procedure* find_procedure(const std::string& id) {
  for (const Procedure* p : registry())
    if (p->id == id) return p;
  return nullptr;
}

std::vector<const Procedure*> catalog_procedures() {
  std::vector<const Procedure*> procs = registry();
  const auto& cats = category_order();
  auto cat_index = [&](const std::string& c) {
    auto it = std::find(cats.begin(), cats.end(), c);
    return it == cats.end() ? cats.size() : static_cast<std::size_t>(it - cats.begin());
  };
  std::stable_sort(procs.begin(), procs.end(), [&](const Procedure* a, const Procedure* b) {
    const auto ca = cat_index(a->category), cb = cat_index(b->category);
    if (ca != cb) return ca < cb;
    if (a->category != b->category) return a->category < b->category;
    if (a->subcategory != b->subcategory) return a->subcategory < b->subcategory;
    if (a->order != b->order) return a->order < b->order;
    return a->name < b->name;
  });
  return procs;
}

const std::vector<std::string>& category_order() {
  static const std::vector<std::string> order = {
      "Means",          "Proportions",         "Correlation", "Regression",          "Survival",
      "Rates and Counts", "Variances",          "Nonparametric", "Group-Sequential",  "Agreement",
      "Diagnostic Accuracy", "Confidence Intervals", "Cluster-Randomized", "Simulation"};
  return order;
}

}  // namespace ssc
