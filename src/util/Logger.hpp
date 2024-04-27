#ifndef OPEN_AICG2_PLUS_UTIL_LOGGER_HPP
#define OPEN_AICG2_PLUS_UTIL_LOGGER_HPP

#include <fmt/core.h>

namespace detail
{
void log_impl(std::string);
} // detail

template<typename ...Ts>
[[noreturn]]
void log_fatal(fmt::format_string<Ts...> fmt, Ts&&... args)
{
    using namespace std::literals::string_literals;
    detail::log_impl("[fatal] "s + fmt::format(std::move(fmt), std::forward<Ts>(args)...));
    std::terminate();
}

template<typename ...Ts>
void log_info(fmt::format_string<Ts...> fmt, Ts&&... args)
{
    using namespace std::literals::string_literals;
    detail::log_impl(fmt::format(std::move(fmt), std::forward<Ts>(args)...));
}

template<typename ...Ts>
void log_assert(bool cond, fmt::format_string<Ts...> fmt, Ts&&... args)
{
    using namespace std::literals::string_literals;
    if( ! cond)
    {
        detail::log_impl("[assert] "s + fmt::format(std::move(fmt), std::forward<Ts>(args)...));
        std::terminate();
    }
    return;
}

#endif// OPEN_AICG2_PLUS_UTIL_LOGGER_HPP
