// Full-text search over the procedure catalog.
//
// Every procedure is indexed by its name, keywords, topic, id, parameter
// labels, summary, definitions and references, each field with its own
// weight.  Text is normalized (case, accents, Greek letters, hyphenated
// compounds, plurals) and phrases from a table of synonyms and abbreviations
// map to shared group terms.  A query word matches index terms exactly, as a
// prefix (a word still being typed) or within a small edit distance (a typo);
// a trigram index over the vocabulary finds the edit-distance candidates
// without scanning it.  Procedures are ranked with BM25 over the weighted
// fields, and every distinctive query word must match unless none does.
#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace ssc::search {

struct Hit {
  std::string id;
  double score = 0;
  std::vector<std::string> terms;  // words of the procedure that matched the query
};

struct Result {
  std::vector<Hit> hits;
  std::string query;      // the query as read (a long text: its first 2,000 bytes)
  std::string corrected;  // the query with misspelled words replaced ("" when unchanged)
  bool partial = false;   // no procedure matches every word: closest procedures instead
};

/// Procedures ranked for a free-text query (at most `limit`; 0 = all matches).
Result query(const std::string& text, std::size_t limit = 20);

/// Normalized index terms of a text: lowercase ASCII with accents removed and
/// Greek letters spelled out, plurals reduced, stop words dropped, and the
/// joined form of hyphenated words added ("log-rank" -> log, rank, logrank).
std::vector<std::string> terms(const std::string& text);

/// Damerau-Levenshtein distance (optimal string alignment), or limit + 1 once it exceeds `limit`.
int edit_distance(const std::string& a, const std::string& b, int limit);

}  // namespace ssc::search
