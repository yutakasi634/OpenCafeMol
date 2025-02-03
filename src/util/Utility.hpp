#ifndef OPEN_AICG2_PLUS_UTILITY_HPP
#define OPEN_AICG2_PLUS_UTILITY_HPP

#include <toml.hpp>

#include <iostream>
#include <string>

namespace Utility
{

// ----------------------------------------------------------------------------
// file io

template<typename T>
void write_as_bytes(std::ostream& os, const T& v) noexcept
{
    using Type = typename std::remove_reference<T>::type;
    os.write(reinterpret_cast<const char*>(std::addressof(v)), sizeof(Type));
    return;
}

std::string get_file_suffix(const std::string& filename);
std::string erase_space(std::string&& str);
void clear_file(const std::string& filename);

// ----------------------------------------------------------------------------
// parse toml file

const toml::value& find_either(
        const toml::value& v, const std::string& key1, const std::string& key2);

// find_parameter is a utility function to support the following functionality.
//
// ```toml
// [[forcefields.global]]
// env.sigma = 1.0
// parameters = [
//     {index = 1, sigma = "sigma"},
//     {index = 2, sigma = 2.0},
//     # ...
// ]
// ```
//
// First, it searches `params` with the `name`. If the corresponding value is a
// string, it searches `env` with the string.
template<typename T>
T find_parameter(const toml::value& params, const toml::value& env,
                 const std::string& name)
{
    static_assert(!std::is_same<T, std::string>::value,
                  "string value cannot be aliased");

    if(!params.is_table() || !params.contains(name))
    {
        const toml::error_info err = toml::make_error_info(
                "value " + name + " does not exists", params, "in this table");
        throw std::out_of_range(toml::format_error(err));
    }
    const toml::value& p = params.at(name);
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_empty())
        {
            const toml::error_info err = toml::make_error_info(
                    "named variable \"" + var + "\" used but no env is defined",
                    params, "used here");
            throw std::out_of_range(toml::format_error(err));
        }
        if(!env.is_table() || !env.contains(var))
        {
            const toml::error_info err = toml::make_error_info(
                    "named variable \"" + var + "\" does not exists",
                    env, "in this table");
            throw std::out_of_range(toml::format_error(err));
        }
        return toml::find<T>(env, var);
    }
    return toml::get<T>(p);
}

template<typename T>
T find_parameter_either(const toml::value& params, const toml::value& env,
                        const std::string& name1,  const std::string& name2)
{
    static_assert(!std::is_same<T, std::string>::value,
                  "string value cannot be aliased");

    if(!params.is_table() || (!params.contains(name1) && !params.contains(name2)))
    {
        const toml::error_info err = toml::make_error_info(
            "value " + name1 + " and " + name2 +
            " does not exists", params, "in this table");
        throw std::out_of_range(toml::format_error(err));
    }
    else if(params.contains(name1) && params.contains(name2))
    {
        const toml::error_info err = toml::make_error_info(
            "key duplicates.",
            params.at(name1), "here",
            params.at(name2), "this conflicts with the above value definition");
        throw std::out_of_range(toml::format_error(err));
    }
    toml::value p;
    if     (params.contains(name1)) { p = params.at(name1); }
    else if(params.contains(name2)) { p = params.at(name2); }
    else
    {
        const toml::error_info err = toml::make_error_info(
            "both keys, \"" + name1 + "\" and \"" + name2 + "\", are not found.",
            params, "in this table");
        throw std::runtime_error(toml::format_error(err));
    }

    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_empty())
        {
            const toml::error_info err = toml::make_error_info(
                "named variable \"" + var + "\" used but no env is defined",
                params, "used here");
            throw std::out_of_range(toml::format_error(err));
        }
        if(!env.is_table() || !env.contains(var))
        {
            const toml::error_info err = toml::make_error_info(
                "named variable \"" + var + "\" does not exists",
                env, "in this table");
            throw std::out_of_range(toml::format_error(err));
        }
        return toml::find<T>(env, var);
    }
    return toml::get<T>(p);
}

// find_parameter with an optional value, opt.
template<typename T>
T find_parameter_or(const toml::value& params, const toml::value& env,
                    const std::string& name, const T& opt) noexcept
{
    static_assert(!std::is_same<T, std::string>::value,
                  "string value cannot be aliased");

    if(!params.is_table() || !params.contains(name))
    {
        return opt;
    }
    const toml::value& p = params.at(name);
    if(p.is_string())
    {
        return toml::find_or(env, p.as_string(), opt);
    }
    return toml::get_or(p, opt);
}

// ----------------------------------------------------------------------------
// handling container

template<typename T>
bool contains(std::vector<std::pair<T, T>> pair_list, std::pair<T, T> query)
{
    for(const auto& pair : pair_list)
    {
        if(pair.first == query.first && pair.second == query.second)
        {
            return true;
        }
        if(pair.second == query.first && pair.first == query.second)
        {
            return true;
        }
    }
    return false;
}

template<std::size_t I, typename Tuple, typename F>
void for_each_impl(Tuple& t, F func)
{
    if constexpr(I == std::tuple_size_v<Tuple>)
    {
        return;
    }
    else
    {
        func(I, std::get<I>(t));
        return for_each_impl<I+1>(t, std::move(func));
    }
}
template<std::size_t I, typename Tuple, typename F>
void for_each_impl(const Tuple& t, F func)
{
    if constexpr(I == std::tuple_size_v<Tuple>)
    {
        return;
    }
    else
    {
        func(I, std::get<I>(t));
        return for_each_impl<I+1>(t, std::move(func));
    }
}
template<typename Tuple, typename F>
void for_each(Tuple& t, F func)
{
    return for_each_impl<0>(t, std::move(func));
}
template<typename Tuple, typename F>
void for_each(const Tuple& t, F func)
{
    return for_each_impl<0>(t, std::move(func));
}

} // namespace Utility

#endif // OPEN_AICG2_PLUS_UTILITY_HPP
