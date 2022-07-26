#ifndef OPEN_AICG2_PLUS_UTILITY_HPP
#define OPEN_AICG2_PLUS_UTILITY_HPP

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
    if (v.contains(key1)) { return v.at(key1); }
    if(v.contains(key2))  { return v.at(key2); }

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

#endif // OPEN_AICG2_PLUS_UTILITY_HPP
