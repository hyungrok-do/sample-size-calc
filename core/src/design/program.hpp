#pragma once
#include "plan.hpp"

namespace ssc::design {
json::Value compile_program(const json::Value &request);
json::Value execute_program(const json::Value &plan, Token *token);
} // namespace ssc::design
