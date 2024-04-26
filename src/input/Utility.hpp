#ifndef OPEN_AICG2_PLUS_INPUT_UTILITY_HPP
#define OPEN_AICG2_PLUS_INPUT_UTILITY_HPP

#include <toml11/toml.hpp>

void merge_toml_tables(toml::value& table, const toml::value& other);
void expand_include(toml::value& v);

void add_offset(std::size_t& index, const toml::value& offset);
void add_offset(std::pair<std::size_t, std::size_t>& indices, const toml::value& offset);

void add_offset(std::array<std::size_t, 3>& indices, const toml::value& offset);
void add_offset(std::array<std::size_t, 4>& indices, const toml::value& offset);

// This check all the keys in a table are found in a list.
//     If there is a key that is not found in the range, it warns about the
// corresponding value will be ignored.
//
// Use it as the following.
// ```cpp
// check_keys_available(table, {"foo", "bar", "baz"});
// ```
bool check_keys_available(const toml::value& table,
                          std::initializer_list<std::string> list);

#endif // OPEN_AICG2_PLUS_INPUT_UTILITY_HPP
