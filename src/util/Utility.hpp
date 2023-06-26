#ifndef OPEN_AICG2_PLUS_UTILITY_HPP
#define OPEN_AICG2_PLUS_UTILITY_HPP

#include <string>
#include <algorithm>

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

std::string get_file_suffix(const std::string& filename)
{
    const std::size_t file_suffix_from = filename.rfind(".");
    if(file_suffix_from == std::string::npos)
    {
        throw std::runtime_error(
                "[error] There is no file extension in " + filename + "."
                " The file type can not be specified.");
    }
    const std::size_t file_suffix_len  = filename.length() - file_suffix_from;
    return filename.substr(file_suffix_from, file_suffix_len);
}

std::string erase_space(std::string&& str)
{
    const auto new_end = std::remove_if(str.begin(), str.end(),
                             [](const char x){ return std::isspace(x); });
    str.erase(new_end, str.end());
    return str;
}

void clear_file(const std::string& filename)
{
    std::ofstream ofs(filename);
    if(not ofs.good())
    {
        throw std::runtime_error("file open error : " + filename);
    }
    ofs.close();
    return;
}

// ----------------------------------------------------------------------------
// parse toml file

const toml::value& find_either(
        const toml::value& v, const std::string& key1, const std::string& key2)
{
    // A functor to find a value that corresponds to either of the key.
    // If both key exists, throw an error.
    if(v.contains(key1) && v.contains(key2) != 0)
    {
        std::cerr << toml::format_error("[error] key duplicates.", v.at(key1), "here", v.at(key2),
                                        "this conflicts with the above value definition")
                 << std::endl;
    }
    if(v.contains(key1)) { return v.at(key1); }
    if(v.contains(key2)) { return v.at(key2); }

    throw std::runtime_error(toml::format_error("both keys, \"" + key1 + "\" and \"" + key2 +
                             "\", are not found.", v, "in this table"));
}

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
        throw std::out_of_range(toml::format_error("[error] value " + name +
            " does not exists", params, "in this table"));
    }
    const toml::value& p = params.at(name);
    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_uninitialized())
        {
            throw std::out_of_range(toml::format_error("[error] named variable \"" +
                var + "\" used but no env is defined", params, "used here"));
        }
        if(!env.is_table() || !env.contains(var))
        {
            throw std::out_of_range(toml::format_error("[error] named variable \"" +
                var + "\" does not exists", env, "in this table"));
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
        throw std::out_of_range(toml::format_error("[error] value " + name1 + " and " + name2 +
            " does not exists", params, "in this table"));
    }
    else if(params.contains(name1) && params.contains(name2))
    {
        std::cerr << toml::format_error("[error] key duplicates.",
                                        params.at(name1), "here", params.at(name2),
                                        "this conflicts with the above value definition")
                 << std::endl;
    }
    toml::value p;
    if     (params.contains(name1)) { p = params.at(name1); }
    else if(params.contains(name2)) { p = params.at(name2); }
    else
    {
        throw std::runtime_error(toml::format_error("both keys, \"" + name1 + "\" and \"" + name2 +
                             "\", are not found.", params, "in this table"));
    }

    if(p.is_string())
    {
        // search inside of `env`
        const std::string& var = p.as_string();
        if(env.is_uninitialized())
        {
            throw std::out_of_range(toml::format_error("[error] named variable \"" +
                var + "\" used but no env is defined", params, "used here"));
        }
        if(!env.is_table() || !env.contains(var))
        {
            throw std::out_of_range(toml::format_error("[error] named variable \"" +
                var + "\" does not exists", env, "in this table"));
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

} // namespace Utility

#endif // OPEN_AICG2_PLUS_UTILITY_HPP
