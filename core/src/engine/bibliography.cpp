#include "engine/bibliography.hpp"

#include <stdexcept>

namespace ssc::bibliography {
namespace {
struct Record {
  std::string key, citation, bibtex, source, export_url;
};
struct Mapping {
  std::string id;
  std::vector<std::string> keys;
};

#include "engine/bibliography_data.inc"

const Record& find(const std::string& key) {
  for (const auto& record : records())
    if (record.key == key) return record;
  throw std::logic_error("Unknown bibliography key: " + key);
}
}  // namespace

void attach(Procedure& procedure) {
  for (const auto& mapping : mappings()) {
    if (mapping.id != procedure.id) continue;
    procedure.reference_keys = mapping.keys;
    procedure.references.clear();
    for (const auto& key : mapping.keys) procedure.references.push_back(find(key).citation);
    return;
  }
}

json::Value to_json(const Procedure& procedure) {
  json::Value out = json::Value::array();
  for (const auto& key : procedure.reference_keys) {
    const auto& record = find(key);
    json::Value item = json::Value::object();
    item.set("key", record.key);
    item.set("citation", record.citation);
    item.set("bibtex", record.bibtex);
    item.set("source", record.source);
    item.set("export_url", record.export_url);
    out.push(item);
  }
  return out;
}
}  // namespace ssc::bibliography
