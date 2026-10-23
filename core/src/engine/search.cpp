#include "engine/search.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <set>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

#include "engine/procedure.hpp"

namespace ssc::search {
namespace {

// ===================================================================== normalization
// Decode the UTF-8 code point at s[i] and advance i (invalid bytes give U+FFFD).
char32_t next_code_point(const std::string& s, std::size_t& i) {
  auto byte = [&](std::size_t k) { return static_cast<unsigned char>(s[i + k]); };
  auto cont = [&](std::size_t k) { return i + k < s.size() && (byte(k) & 0xC0) == 0x80; };
  const unsigned char c = byte(0);
  if (c < 0x80) {
    ++i;
    return c;
  }
  if ((c & 0xE0) == 0xC0 && cont(1)) {
    const char32_t cp = (char32_t(c & 0x1F) << 6) | (byte(1) & 0x3F);
    i += 2;
    return cp;
  }
  if ((c & 0xF0) == 0xE0 && cont(1) && cont(2)) {
    const char32_t cp = (char32_t(c & 0x0F) << 12) | (char32_t(byte(1) & 0x3F) << 6) | (byte(2) & 0x3F);
    i += 3;
    return cp;
  }
  if ((c & 0xF8) == 0xF0 && cont(1) && cont(2) && cont(3)) {
    const char32_t cp = (char32_t(c & 0x07) << 18) | (char32_t(byte(1) & 0x3F) << 12) |
                        (char32_t(byte(2) & 0x3F) << 6) | (byte(3) & 0x3F);
    i += 4;
    return cp;
  }
  ++i;
  return 0xFFFD;
}

// ASCII spelling of a non-ASCII letter or digit ("" for anything else).
const char* spell(char32_t cp) {
  if ((cp >= 0xC0 && cp <= 0xC5) || (cp >= 0xE0 && cp <= 0xE5)) return "a";
  if ((cp >= 0xC8 && cp <= 0xCB) || (cp >= 0xE8 && cp <= 0xEB)) return "e";
  if ((cp >= 0xCC && cp <= 0xCF) || (cp >= 0xEC && cp <= 0xEF)) return "i";
  if ((cp >= 0xD2 && cp <= 0xD6) || cp == 0xD8 || (cp >= 0xF2 && cp <= 0xF6) || cp == 0xF8) return "o";
  if ((cp >= 0xD9 && cp <= 0xDC) || (cp >= 0xF9 && cp <= 0xFC)) return "u";
  switch (cp) {
    case 0xC6: case 0xE6: return "ae";
    case 0xC7: case 0xE7: case 0x0107: case 0x010D: return "c";
    case 0xD0: case 0xF0: return "d";
    case 0xD1: case 0xF1: case 0x0144: case 0x0148: return "n";
    case 0xDD: case 0xFD: case 0xFF: return "y";
    case 0xDE: case 0xFE: return "th";
    case 0xDF: return "ss";
    case 0xD7: return "x";  // multiplication sign: 2×2 -> 2x2
    case 0x0142: case 0x0141: return "l";
    case 0x0161: case 0x0160: case 0x015B: case 0x015F: return "s";
    case 0x017E: case 0x017D: case 0x017A: case 0x017C: return "z";
    case 0x0151: return "o";
    case 0x0171: return "u";
    case 0x0119: return "e";
    case 0x0105: return "a";
    case 0x0159: return "r";
    case 0x011F: return "g";
    case 0x0131: return "i";
    case 0x0152: case 0x0153: return "oe";
    // superscript and subscript digits
    case 0xB9: case 0x2081: return "1";
    case 0xB2: case 0x2082: return "2";
    case 0xB3: case 0x2083: return "3";
    case 0x2070: case 0x2080: return "0";
    // Greek letters, spelled out (σ and µ have two code points each)
    case 0x03B1: case 0x0391: return "alpha";
    case 0x03B2: case 0x0392: return "beta";
    case 0x03B3: case 0x0393: return "gamma";
    case 0x03B4: case 0x0394: return "delta";
    case 0x03B5: case 0x0395: return "epsilon";
    case 0x03B6: case 0x0396: return "zeta";
    case 0x03B7: case 0x0397: return "eta";
    case 0x03B8: case 0x0398: return "theta";
    case 0x03B9: case 0x0399: return "iota";
    case 0x03BA: case 0x039A: return "kappa";
    case 0x03BB: case 0x039B: return "lambda";
    case 0x03BC: case 0x039C: case 0xB5: return "mu";
    case 0x03BD: case 0x039D: return "nu";
    case 0x03BE: case 0x039E: return "xi";
    case 0x03C0: case 0x03A0: return "pi";
    case 0x03C1: case 0x03A1: return "rho";
    case 0x03C2: case 0x03C3: case 0x03A3: return "sigma";
    case 0x03C4: case 0x03A4: return "tau";
    case 0x03C6: case 0x03A6: case 0x03D5: return "phi";
    case 0x03C7: case 0x03A7: return "chi";
    case 0x03C8: case 0x03A8: return "psi";
    case 0x03C9: case 0x03A9: return "omega";
    default: break;
  }
  static const char* const digits[] = {"4", "5", "6", "7", "8", "9"};
  if (cp >= 0x2074 && cp <= 0x2079) return digits[cp - 0x2074];
  if (cp >= 0x2084 && cp <= 0x2089) return digits[cp - 0x2084];
  return "";
}

bool is_hyphen(char32_t cp) { return cp == '-' || (cp >= 0x2010 && cp <= 0x2013) || cp == 0x2212; }
bool is_apostrophe(char32_t cp) { return cp == '\'' || cp == 0x2019 || cp == 0x02BC; }

struct Word {
  std::string text;     // lowercase ASCII
  bool joined = false;  // a hyphen joins it to the previous word
};

std::vector<Word> split_words(const std::string& s) {
  std::vector<Word> out;
  Word cur;
  bool hyphen_pending = false;
  auto flush = [&]() {
    if (!cur.text.empty()) out.push_back(cur);
    cur = Word();
  };
  auto append = [&](const char* text) {
    if (cur.text.empty()) cur.joined = hyphen_pending;
    hyphen_pending = false;
    cur.text += text;
  };
  for (std::size_t i = 0; i < s.size();) {
    const char32_t cp = next_code_point(s, i);
    if (cp < 0x80 && std::isalnum(static_cast<int>(cp))) {
      const char c[2] = {static_cast<char>(std::tolower(static_cast<int>(cp))), 0};
      append(c);
    } else if (const char* sp = spell(cp); *sp) {
      append(sp);
    } else if (is_apostrophe(cp)) {
      continue;  // Cohen's -> cohens
    } else if (is_hyphen(cp) && !cur.text.empty()) {
      flush();
      hyphen_pending = true;
    } else {
      flush();
      hyphen_pending = false;
    }
  }
  flush();
  return out;
}

const std::unordered_set<std::string>& stop_words() {
  static const std::unordered_set<std::string> words = {
      "a",     "an",   "the",     "of",   "for",   "and",   "or",     "in",   "on",   "at",    "to",
      "with",  "by",   "from",    "is",   "are",   "be",    "as",     "vs",   "versus", "via", "using",
      "use",   "between", "under", "into", "its",  "it",    "this",   "that", "these", "those", "what",
      "which", "how",  "when",    "i",    "we",    "my",    "our",    "do",   "does",  "need",  "want",
      "calculate", "calculation", "compute", "find", "get", "per", "than", "then", "if", "can", "will",
      "compare", "compared", "comparing", "comparison"};
  return words;
}

bool ends_with(const std::string& s, const char* tail) {
  const std::size_t n = std::char_traits<char>::length(tail);
  return s.size() >= n && s.compare(s.size() - n, n, tail) == 0;
}

// Plural and similar endings (the same reduction is applied to the index and the queries).
std::string stem(const std::string& w) {
  const std::size_t n = w.size();
  if (n <= 3 || !std::isalpha(static_cast<unsigned char>(w.back()))) return w;
  if (ends_with(w, "yses")) return w.substr(0, n - 2) + "is";  // analyses -> analysis
  if (ends_with(w, "ies") && n > 4) return w.substr(0, n - 3) + "y";
  if (ends_with(w, "sses")) return w.substr(0, n - 2);
  if (ends_with(w, "ches") || ends_with(w, "shes") || ends_with(w, "xes")) return w.substr(0, n - 2);
  if (ends_with(w, "ss") || ends_with(w, "us") || ends_with(w, "is") || ends_with(w, "ys")) return w;
  if (w.back() == 's') return w.substr(0, n - 1);
  return w;
}

struct Token {
  std::string term;     // index term
  std::string surface;  // readable form (the word, or the hyphenated compound)
  bool compound = false;
};

// Terms of a text, with the joined form of each run of hyphenated words.
std::vector<Token> analyze(const std::string& text) {
  const std::vector<Word> ws = split_words(text);
  std::vector<Token> out;
  for (std::size_t i = 0; i < ws.size();) {
    std::size_t j = i + 1;
    while (j < ws.size() && ws[j].joined) ++j;
    for (std::size_t k = i; k < j; ++k)
      if (!stop_words().count(ws[k].text)) out.push_back({stem(ws[k].text), ws[k].text, false});
    if (j - i > 1) {
      std::string joined, surface;
      for (std::size_t k = i; k < j; ++k) {
        joined += ws[k].text;
        surface += (k > i ? "-" : "") + ws[k].text;
      }
      out.push_back({stem(joined), surface, true});
      for (std::size_t k = i; j - i > 2 && k + 1 < j; ++k)
        out.push_back({stem(ws[k].text + ws[k + 1].text), ws[k].text + "-" + ws[k + 1].text, true});
    }
    i = j;
  }
  return out;
}

// ===================================================================== synonyms
// Phrases that mean the same thing for searching; each group becomes one shared
// term ("~" + its first phrase) wherever any of its phrases occurs.
const char* const SYNONYMS[] = {
    "noninferiority|non inferiority|ni",
    "equivalence|tost|two one sided test|two one sided tests",
    "bioequivalence|bioavailability|average bioequivalence",
    "confidence interval|ci|precision|margin of error|interval width",
    "anova|analysis of variance",
    "ancova|analysis of covariance|covariate adjusted|covariate adjustment|adjusted for covariates",
    "auc|roc|roc curve|area under the curve|area under curve|receiver operating characteristic",
    "icc|intraclass correlation|intra class correlation",
    "relative risk|risk ratio|rr",
    "hazard ratio|hr",
    "kaplan meier|km",
    "group sequential|gsd|interim analysis|interim analyses|interim look|interim monitoring|early stopping",
    "ttest|t test|student t|students t",
    "ztest|z test",
    "chi square|chisquare|chi squared|chisq|chi2",
    "mann whitney|rank sum|wilcoxon rank sum|u test|wilcoxon mann whitney",
    "signed rank|wilcoxon signed rank",
    "crossover|cross over",
    "randomized controlled trial|rct|randomized trial|parallel group",
    "survival|time to event|event time|failure time",
    "logrank|log rank",
    "poisson|count data|event count|incidence rate",
    "negative binomial|overdispersion|overdispersed|recurrent event|exacerbation",
    "agreement|reliability|concordance|inter rater|interrater|inter observer|interobserver|rater",
    "diagnostic|screening|diagnostic accuracy|diagnostic test",
    "longitudinal|repeated measure|repeated measurement|within subject|growth curve|follow up visit",
    "dose response|trend|ordered group|ordinal group",
    "phase ii|phase 2|single arm",
    "cluster randomized|cluster trial|group randomized|community trial|cluster",
    "proportion|binomial|percentage|binary outcome|dichotomous|response rate",
    "variance|standard deviation|variability|dispersion",
    "paired|matched pair|matched|before after|pre post",
    "one sample|single sample|single group|one group",
    "two sample|two group|independent sample|independent group|two independent",
    "correlation|pearson|association",
    "effect size|cohen d|cohens d|standardized difference|detectable difference",
    "sensitivity|specificity|true positive rate|true negative rate",
    "stepped wedge|staggered rollout",
    "kappa|cohen kappa|cohens kappa",
    "fisher exact|fishers exact|fisher test",
    "mcnemar|mc nemar",
    "dropout|attrition|loss to follow up|lost to follow up",
    "logistic|logit",
    "slope|rate of change|change over time",
};

struct Synonyms {
  std::vector<std::pair<std::vector<std::string>, std::string>> phrases;  // longest first
  Synonyms() {
    for (const char* group : SYNONYMS) {
      std::string text(group), name;
      std::size_t start = 0;
      while (start <= text.size()) {
        const std::size_t bar = std::min(text.find('|', start), text.size());
        std::vector<std::string> phrase;
        for (const Token& t : analyze(text.substr(start, bar - start)))
          if (!t.compound) phrase.push_back(t.term);
        if (!phrase.empty()) {
          if (name.empty()) {
            name = "~";
            for (const auto& p : phrase) name += (name.size() > 1 ? "_" : "") + p;
          }
          phrases.emplace_back(std::move(phrase), name);
        }
        start = bar + 1;
      }
    }
    std::stable_sort(phrases.begin(), phrases.end(),
                     [](const auto& a, const auto& b) { return a.first.size() > b.first.size(); });
  }
  // The group and phrase length that starts at terms[i] (longest phrase), if any.
  std::pair<const std::string*, std::size_t> at(const std::vector<std::string>& terms, std::size_t i) const {
    for (const auto& [phrase, name] : phrases) {
      if (i + phrase.size() > terms.size()) continue;
      if (std::equal(phrase.begin(), phrase.end(), terms.begin() + static_cast<std::ptrdiff_t>(i)))
        return {&name, phrase.size()};
    }
    return {nullptr, 0};
  }
};

const Synonyms& synonyms() {
  static const Synonyms s;
  return s;
}

// ===================================================================== the index
enum Field { F_NAME, F_KEYWORDS, F_TOPIC, F_ID, F_PARAMS, F_SUMMARY, F_DEFINITIONS, F_REFERENCES, F_COUNT };
constexpr double FIELD_WEIGHT[F_COUNT] = {3.0, 2.5, 2.0, 2.0, 0.7, 1.0, 0.5, 0.4};
constexpr double K1 = 1.2;           // BM25 term-frequency saturation
constexpr double COMMON_SHARE = 0.6;  // words in more procedures than this are optional in queries

struct Posting {
  std::uint32_t doc;
  double weight;     // field-weighted term frequency
  bool cited_only;   // only in the references: adds to the score but does not satisfy a query word
};

struct Index {
  std::vector<std::string> ids, names;
  std::vector<std::string> name_text;  // names as phrase_of() reads them, for the phrase bonus
  std::vector<std::string> vocab;      // sorted
  std::unordered_map<std::string, std::uint32_t> lookup;
  std::vector<std::vector<Posting>> postings;
  std::vector<double> idf;
  std::vector<std::string> surface;
  std::unordered_map<std::string, std::vector<std::uint32_t>> grams;  // trigram -> terms
};

std::vector<std::string> trigrams(const std::string& term) {
  const std::string padded = "$" + term + "$";
  std::vector<std::string> out;
  for (std::size_t i = 0; i + 3 <= padded.size(); ++i) out.push_back(padded.substr(i, 3));
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// The words of a text for comparing phrases: synonym phrases become their
// group term and the generic word "test" is left out, so that the name
// "Two-Sample T-Tests" and the query "two sample t test" read alike.
std::string phrase_of(const std::string& text) {
  std::vector<std::string> plain;
  for (const Token& t : analyze(text))
    if (!t.compound) plain.push_back(t.term);
  std::string out;
  for (std::size_t i = 0; i < plain.size();) {
    const auto [group, len] = synonyms().at(plain, i);
    const std::string& word = group ? *group : plain[i];
    if (word != "test") out += (out.empty() ? "" : " ") + word;
    i += group ? len : 1;
  }
  return out;
}

void build(Index& idx) {
  const std::vector<const Procedure*> procs = catalog_procedures();

  std::map<std::string, std::map<std::uint32_t, double>> tf;  // term -> doc -> weighted frequency
  std::map<std::string, std::set<std::uint32_t>> described;   // term -> docs using it outside the references
  std::map<std::string, std::map<std::string, int>> forms;     // term -> surface -> count
  for (std::uint32_t d = 0; d < procs.size(); ++d) {
    const Procedure& p = *procs[d];
    idx.ids.push_back(p.id);
    idx.names.push_back(p.name);
    std::string fields[F_COUNT];
    fields[F_NAME] = p.name;
    for (const auto& k : p.keywords) fields[F_KEYWORDS] += k + " ; ";
    fields[F_TOPIC] = p.category + " ; " + p.subcategory;
    fields[F_ID] = p.id;
    std::replace(fields[F_ID].begin(), fields[F_ID].end(), '.', ' ');
    for (const auto& s : p.solve_for) fields[F_PARAMS] += s.label + " ; ";
    for (const auto& q : p.params) {
      fields[F_PARAMS] += q.label + " ; " + q.symbol + " ; ";
      for (const auto& o : q.options) fields[F_PARAMS] += o.label + " ; ";
    }
    fields[F_SUMMARY] = p.summary;
    for (const auto& [sym, text] : p.definitions) fields[F_DEFINITIONS] += sym + " " + text + " ; ";
    for (const auto& r : p.references) fields[F_REFERENCES] += r + " ; ";
    idx.name_text.push_back(" " + phrase_of(p.name) + " ");
    for (int f = 0; f < F_COUNT; ++f) {
      const std::vector<Token> tokens = analyze(fields[f]);
      std::map<std::string, int> count;
      std::vector<std::string> plain;
      for (const Token& t : tokens) {
        ++count[t.term];
        ++forms[t.term][t.surface];
        if (!t.compound) plain.push_back(t.term);
      }
      std::unordered_set<std::string> groups;
      for (std::size_t i = 0; i < plain.size(); ++i)
        if (const auto [group, len] = synonyms().at(plain, i); group) groups.insert(*group);
      for (const auto& g : groups) count[g] = std::max(count[g], 1);
      for (const auto& [term, c] : count) {
        tf[term][d] += FIELD_WEIGHT[f] * std::min(c, 3);
        if (f != F_REFERENCES) described[term].insert(d);
      }
    }
  }
  const double n_docs = static_cast<double>(procs.size());
  for (const auto& [term, docs] : tf) {
    const auto id = static_cast<std::uint32_t>(idx.vocab.size());
    idx.vocab.push_back(term);
    idx.lookup[term] = id;
    std::vector<Posting> list;
    for (const auto& [doc, w] : docs) list.push_back({doc, w, described[term].count(doc) == 0});
    idx.postings.push_back(std::move(list));
    const double df = static_cast<double>(docs.size());
    idx.idf.push_back(std::log(1 + (n_docs - df + 0.5) / (df + 0.5)));
    std::string best = term;
    int most = 0;
    for (const auto& [form, c] : forms[term])
      if (c > most) {
        most = c;
        best = form;
      }
    idx.surface.push_back(term[0] == '~' ? term.substr(1) : best);
    if (term[0] != '~' && term.size() >= 3)
      for (const auto& g : trigrams(term)) idx.grams[g].push_back(id);
  }
}

const Index& the_index() {
  static Index idx;
  static std::once_flag once;
  std::call_once(once, [] { build(idx); });
  return idx;
}

// ===================================================================== queries
struct Alternative {
  std::uint32_t term;
  double weight;
};

struct Clause {
  std::vector<Alternative> alts;
  bool required = true;
  std::string text;       // the query words this clause stands for
  std::string corrected;  // their spelling correction ("" if none)
};

// Terms within edit distance of `word` (typos), found through shared trigrams.
std::vector<std::pair<std::uint32_t, int>> fuzzy(const Index& idx, const std::string& word) {
  const int k = word.size() >= 8 ? 2 : 1;
  const std::vector<std::string> grams = trigrams(word);
  std::unordered_map<std::uint32_t, int> shared;
  for (const auto& g : grams)
    if (const auto it = idx.grams.find(g); it != idx.grams.end())
      for (std::uint32_t t : it->second) ++shared[t];
  const int need = std::max(1, static_cast<int>(grams.size()) - 4 * k);
  std::vector<std::pair<std::uint32_t, int>> out;
  for (const auto& [t, c] : shared) {
    if (c < need) continue;
    const std::string& term = idx.vocab[t];
    if (std::abs(static_cast<int>(term.size()) - static_cast<int>(word.size())) > k) continue;
    const int d = edit_distance(word, term, k);
    if (d >= 1 && d <= k) out.emplace_back(t, d);
  }
  std::sort(out.begin(), out.end(), [&](const auto& a, const auto& b) {
    if (a.second != b.second) return a.second < b.second;
    if (idx.postings[a.first].size() != idx.postings[b.first].size())
      return idx.postings[a.first].size() > idx.postings[b.first].size();
    return a.first < b.first;
  });
  return out;
}

// `typing`: the last word of a query that does not end in a space, which may
// still be incomplete.
Clause word_clause(const Index& idx, const std::string& term, const std::string& text, bool typing) {
  Clause c;
  c.text = text;
  const auto exact = idx.lookup.find(term);
  const bool known = exact != idx.lookup.end();
  if (known) c.alts.push_back({exact->second, 1.0});
  bool prefixed = false;
  // Prefixes: a word that is not in the index is taken as unfinished.  A known
  // word only while it is being typed ("log" on the way to "logrank"), and then
  // its longer matches never count more than the word itself: rare terms such
  // as "two-stage" or "testing" would otherwise outrank "two" or "test".
  if (term.size() >= 3 && (!known || typing)) {
    auto it = std::lower_bound(idx.vocab.begin(), idx.vocab.end(), term);
    for (int n = 0; it != idx.vocab.end() && it->compare(0, term.size(), term) == 0 && n < 40; ++it, ++n) {
      if (*it == term) continue;
      const auto t = static_cast<std::uint32_t>(it - idx.vocab.begin());
      c.alts.push_back({t, known ? 0.85 * std::min(1.0, idx.idf[exact->second] / idx.idf[t]) : 0.85});
      prefixed = true;
    }
  }
  if (term.size() >= 4 && !known) {  // a typo
    const auto near = fuzzy(idx, term);
    for (const auto& [t, d] : near) c.alts.push_back({t, d == 1 ? 0.75 : 0.55});
    if (!near.empty() && !prefixed) c.corrected = idx.surface[near.front().first];
  }
  return c;
}

std::vector<Clause> parse_query(const Index& idx, const std::string& text) {
  const bool typing = !text.empty() && std::isalnum(static_cast<unsigned char>(text.back()));
  // words, with hyphenated runs replaced by their joined form when the index knows it
  struct Unit {
    std::string term, text;
  };
  std::vector<Unit> units;
  const std::vector<Word> ws = split_words(text);
  for (std::size_t i = 0; i < ws.size();) {
    std::size_t j = i + 1;
    while (j < ws.size() && ws[j].joined) ++j;
    std::string joined, shown;
    for (std::size_t k = i; k < j; ++k) {
      joined += ws[k].text;
      shown += (k > i ? "-" : "") + ws[k].text;
    }
    if (j - i > 1 && idx.lookup.count(stem(joined))) {
      units.push_back({stem(joined), shown});
    } else {
      for (std::size_t k = i; k < j; ++k)
        if (!stop_words().count(ws[k].text)) units.push_back({stem(ws[k].text), ws[k].text});
    }
    i = j;
  }
  // synonym phrases become one clause on their group term; other words one clause each
  std::vector<std::string> seq;
  for (const auto& u : units) seq.push_back(u.term);
  std::vector<Clause> clauses;
  std::vector<std::size_t> first_unit;  // the unit each clause starts at
  for (std::size_t i = 0; i < units.size();) {
    const auto [group, len] = synonyms().at(seq, i);
    first_unit.push_back(i);
    if (group) {
      Clause c;
      for (std::size_t k = i; k < i + len; ++k) c.text += (k > i ? " " : "") + units[k].text;
      if (const auto it = idx.lookup.find(*group); it != idx.lookup.end()) c.alts.push_back({it->second, 1.0});
      if (len == 1)  // the word itself counts a little more than its synonyms
        if (const auto it = idx.lookup.find(units[i].term); it != idx.lookup.end())
          c.alts.push_back({it->second, 1.15});
      clauses.push_back(std::move(c));
      i += len;
    } else {
      clauses.push_back(word_clause(idx, units[i].term, units[i].text, typing && i + 1 == units.size()));
      ++i;
    }
  }
  // two neighboring words that the index also knows joined ("log rank", "two sample") are
  // each satisfied by the joined term too
  for (std::size_t c = 0; c + 1 < clauses.size(); ++c) {
    const Unit& a = units[first_unit[c]];
    const Unit& b = units[first_unit[c + 1]];
    if (const auto it = idx.lookup.find(stem(a.text + b.text)); it != idx.lookup.end()) {
      clauses[c].alts.push_back({it->second, 1.0});
      clauses[c + 1].alts.push_back({it->second, 1.0});
    }
  }
  const double n_docs = static_cast<double>(idx.ids.size());
  for (Clause& c : clauses) {
    bool common = !c.alts.empty();
    for (const auto& a : c.alts) common = common && idx.postings[a.term].size() > COMMON_SHARE * n_docs;
    c.required = !common;
  }
  return clauses;
}

}  // namespace

// ===================================================================== public API
std::vector<std::string> terms(const std::string& text) {
  std::vector<std::string> out;
  for (const Token& t : analyze(text)) out.push_back(t.term);
  return out;
}

int edit_distance(const std::string& a, const std::string& b, int limit) {
  const std::size_t n = a.size(), m = b.size();
  if (static_cast<int>(n > m ? n - m : m - n) > limit) return limit + 1;
  std::vector<int> prev2(m + 1), prev(m + 1), cur(m + 1);
  for (std::size_t j = 0; j <= m; ++j) prev[j] = static_cast<int>(j);
  for (std::size_t i = 1; i <= n; ++i) {
    cur[0] = static_cast<int>(i);
    int best = cur[0];
    for (std::size_t j = 1; j <= m; ++j) {
      const int cost = a[i - 1] == b[j - 1] ? 0 : 1;
      cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost});
      if (i > 1 && j > 1 && a[i - 1] == b[j - 2] && a[i - 2] == b[j - 1]) cur[j] = std::min(cur[j], prev2[j - 2] + 1);
      best = std::min(best, cur[j]);
    }
    if (best > limit) return limit + 1;
    std::swap(prev2, prev);
    std::swap(prev, cur);
  }
  return std::min(prev[m], limit + 1);
}

Result query(const std::string& full, std::size_t limit) {
  const Index& idx = the_index();
  // a study is described in a few sentences: of a longer text (a pasted document) the first 2,000 bytes are read
  // (in whole characters) and at most 64 of its words count, so a megabyte cost seconds and gigabytes no more
  std::size_t cut = std::min<std::size_t>(full.size(), 2000);
  while (cut < full.size() && cut > 0 && (static_cast<unsigned char>(full[cut]) & 0xC0) == 0x80) --cut;
  const std::string text = full.substr(0, cut);
  std::vector<Clause> clauses = parse_query(idx, text);
  if (clauses.size() > 64) clauses.resize(64);
  Result res;
  res.query = text;
  if (clauses.empty()) return res;
  const std::size_t n = idx.ids.size();
  std::vector<double> total(n, 0.0);
  std::vector<int> satisfied(n, 0), required_met(n, 0);
  std::vector<std::vector<std::string>> matched(n);
  int required = 0;
  for (const Clause& c : clauses) {
    required += c.required;
    std::vector<double> best(n, 0.0);
    std::vector<std::uint32_t> via(n, 0);
    std::vector<char> described(n, 0);  // matched outside the references (a citation alone does not count)
    for (const auto& a : c.alts)
      for (const Posting& p : idx.postings[a.term]) {
        const double s = a.weight * idx.idf[a.term] * p.weight * (K1 + 1) / (p.weight + K1);
        if (s > best[p.doc]) {
          best[p.doc] = s;
          via[p.doc] = a.term;
        }
        described[p.doc] |= !p.cited_only;
      }
    for (std::size_t d = 0; d < n; ++d) {
      if (best[d] <= 0) continue;
      total[d] += best[d];
      if (std::find(matched[d].begin(), matched[d].end(), idx.surface[via[d]]) == matched[d].end())
        matched[d].push_back(idx.surface[via[d]]);  // each term once
      if (!described[d]) continue;
      ++satisfied[d];
      required_met[d] += c.required;
    }
  }
  // the query with misspelled words replaced
  bool changed = false;
  std::string corrected;
  for (const Clause& c : clauses) {
    corrected += (corrected.empty() ? "" : " ") + (c.corrected.empty() ? c.text : c.corrected);
    changed = changed || !c.corrected.empty();
  }
  if (changed) res.corrected = corrected;

  // a bonus when the (corrected) query reads as a phrase of the name
  const std::string phrase = phrase_of(changed ? corrected : text);
  if (!phrase.empty()) {
    const std::string padded = " " + phrase + " ";
    for (std::size_t d = 0; d < n; ++d) {
      const std::size_t at = idx.name_text[d].find(padded);
      if (at != std::string::npos) total[d] += at == 0 ? 3.0 : 2.0;  // more when the name starts with it
    }
  }

  // a procedure id typed exactly ("prop.two") comes first
  std::string typed;
  for (const char ch : text)
    if (!std::isspace(static_cast<unsigned char>(ch))) typed += static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  for (std::size_t d = 0; d < n; ++d)
    if (idx.ids[d] == typed) total[d] += 100.0;

  std::vector<std::size_t> order;
  for (std::size_t d = 0; d < n; ++d)
    if (satisfied[d] > 0 && required_met[d] == required) order.push_back(d);
  std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return total[a] > total[b]; });
  // Procedures matching only some of the words (most words first) follow: in a ranked list
  // (limit > 0) they fill the room left, so that one procedure matching every word does not hide
  // the others; as a filter (limit 0) only the closest ones, and only when nothing matches fully.
  res.partial = order.empty();
  if (order.empty() || (limit > 0 && order.size() < limit)) {
    std::vector<std::size_t> rest;
    for (std::size_t d = 0; d < n; ++d)
      if (satisfied[d] > 0 && required_met[d] < required) rest.push_back(d);
    std::stable_sort(rest.begin(), rest.end(), [&](std::size_t a, std::size_t b) {
      return satisfied[a] != satisfied[b] ? satisfied[a] > satisfied[b] : total[a] > total[b];
    });
    if (limit == 0 && !rest.empty()) {
      const int most = satisfied[rest[0]];
      rest.erase(std::find_if(rest.begin(), rest.end(), [&](std::size_t d) { return satisfied[d] < most; }), rest.end());
    }
    order.insert(order.end(), rest.begin(), rest.end());
  }
  res.partial = res.partial && !order.empty();
  if (limit > 0 && order.size() > limit) order.resize(limit);
  for (std::size_t d : order) res.hits.push_back({idx.ids[d], total[d], matched[d]});
  return res;
}

}  // namespace ssc::search
