// Parsing of multi-valued numeric inputs.
//
//   "0.8"                 single value
//   "0.8 0.9" / "0.8, 0.9" list (blanks and/or commas)
//   "10 to 100 by 10"     arithmetic range (also "10 to 15", step 1)
//   "5 10 to 30 by 10 50" lists and ranges can be mixed
// Vector inputs (e.g. group means) use the same syntax for one vector; several
// alternative vectors are separated by ';' or '|'.
#pragma once

#include <string>
#include <vector>

namespace ssc::parse {

struct Numbers {
  std::vector<double> values;
  std::string error;  // empty on success
};

struct Vectors {
  std::vector<std::vector<double>> vectors;
  std::string error;
};

/// The longest input text read (10,000 values need far less).
constexpr std::size_t MAX_TEXT = 1000000;

/// Text from a request quoted in a message: at most 40 bytes, cut where a character begins (a field can hold
/// 1,000,000 characters, and a page lays out an unbroken word slowly), with the characters that would act on the
/// message (controls, bidirectional and joining format characters, combining marks) written as \u{...}.
std::string quoted(const std::string& text);

/// The number of characters (code points) of UTF-8 text.
std::size_t characters(const std::string& text);

Numbers numbers(const std::string& text, std::size_t max_values = 10000);
Vectors vectors(const std::string& text, std::size_t max_vectors = 1000);

}  // namespace ssc::parse
