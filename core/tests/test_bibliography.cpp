#include <set>

#include "engine/bibliography.hpp"
#include "engine/runner.hpp"
#include "testing.hpp"

using namespace ssc;

TEST_CASE("bibliography: every registered procedure resolves retrieved records") {
  for (const auto* procedure : ssc::registry()) {
    const auto schema = procedure_json(*procedure);
    const auto& records = schema.get("bibliography");
    CHECK(records.size() > 0);
    CHECK_EQ(records.size(), procedure->references.size());
    std::set<std::string> keys;
    for (std::size_t i = 0; i < records.size(); ++i) {
      const auto& record = records.elements()[i];
      const auto key = record.get("key").as_string();
      CHECK(keys.insert(key).second);
      CHECK_EQ(record.get("citation").as_string(), procedure->references[i]);
      CHECK(record.get("bibtex").as_string().find("{" + key + ",") != std::string::npos);
      CHECK(!record.get("source").as_string().empty());
      CHECK(record.get("export_url").as_string().find("https://") == 0);
    }
  }
}

TEST_CASE("bibliography: calculation results carry the same records as the schema") {
  Token token;
  const auto result = run_request(json::parse(
      R"({"procedure":"means.two.t","solve_for":"n","inputs":{"power":"0.8","mu1":"1","mu2":"0","sigma":"1"}})"), &token);
  CHECK(result.get("ok").as_bool());
  const auto* procedure = find_procedure("means.two.t");
  const auto expected = bibliography::to_json(*procedure);
  CHECK_EQ(result.get("bibliography").size(), expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK_EQ(result.get("bibliography").elements()[i].get("bibtex").as_string(),
             expected.elements()[i].get("bibtex").as_string());
  }
}
