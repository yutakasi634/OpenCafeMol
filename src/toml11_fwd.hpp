#ifndef TOML11_FORWARD_DECLARATION_HPP
#define TOML11_FORWARD_DECLARATION_HPP

#ifdef TOML11_PRESERVE_COMMENTS_BY_DEFAULT
#  define TOML11_DEFAULT_COMMENT_STRATEGY ::toml::preserve_comments
#else
#  define TOML11_DEFAULT_COMMENT_STRATEGY ::toml::discard_comments
#endif

#include <unordered_map>
#include <vector>

namespace toml
{
struct discard_comments;
struct preserve_comments;

template<typename Comment, // discard/preserve_comment
         template<typename ...> class Table,
         template<typename ...> class Array>
class basic_value;
using value = basic_value<TOML11_DEFAULT_COMMENT_STRATEGY, std::unordered_map, std::vector>;

} // toml
#endif//TOML11_FORWARD_DECLARATION_HPP
