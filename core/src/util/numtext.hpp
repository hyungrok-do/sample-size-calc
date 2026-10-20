// Locale-independent conversion between doubles and text.
//
// strtod/printf honour the process locale (LC_NUMERIC), which GUI toolkits
// may set to e.g. de_DE where the decimal separator is ','.  Everything here
// uses the classic "C" locale explicitly, so results never depend on the host.
#pragma once

#include <string>

namespace ssc::text {

/// Parse a floating-point number at s[pos...]; on success sets `out`, advances
/// `pos` past the number and returns true.  Accepts inf/infinity/nan words.
bool parse_double(const std::string& s, std::size_t& pos, double& out);

/// Parse an entire string as one number (surrounding blanks allowed).
bool parse_double(const std::string& s, double& out);

/// printf("%.*g") equivalent in the classic locale.
std::string format_g(double d, int precision);

/// Fixed notation with `decimals` digits after the point (classic locale).
std::string format_fixed(double d, int decimals);

/// Shortest representation that round-trips (≤ 17 significant digits).
std::string format_roundtrip(double d);

}  // namespace ssc::text
