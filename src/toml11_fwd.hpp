#ifndef TOML11_FORWARD_DECLARATION_HPP
#define TOML11_FORWARD_DECLARATION_HPP

#include <unordered_map>
#include <vector>
#include <string>

namespace toml
{

struct type_config;
template<typename type_config>
class basic_value;
using value = basic_value<type_config>;

} // toml
#endif//TOML11_FORWARD_DECLARATION_HPP
