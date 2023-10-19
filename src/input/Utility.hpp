#ifndef OPEN_AICG2_PLUS_INPUT_UTILITY_HPP
#define OPEN_AICG2_PLUS_INPUT_UTILITY_HPP

inline void merge_toml_tables(toml::value& table, const toml::value& other)
{
    assert(table.is_table());
    assert(other.is_table());

    for(const auto& kv : other.as_table())
    {
        if(table.contains(kv.first))
        {
            if(table.at(kv.first).is_table())
            {
                merge_toml_tables(table.at(kv.first), kv.second);
            }
            else
            {
                throw std::runtime_error(toml::format_error("value \"" +
                    kv.first + "\" duplicates.", table, "first defined here",
                    kv.second, "and also defined here"));
            }
        }
        else
        {
            table.as_table().emplace(kv.first, kv.second);
        }
    }
}

inline void expand_include(toml::value& v)
{
    if(!v.is_table() && !v.is_array()) {return;}

    if(v.is_table())
    {
        for(auto& kv : v.as_table())
        {
            if(kv.first == "include") {continue;}
            expand_include(kv.second);
        }

        // expand include in this table
        if(v.contains("include"))
        {
            if(v.at("include").is_array())
            {
                for(auto fname : toml::find<std::vector<std::string>>(v, "include"))
                {
                    std::cerr << "expanding file " << fname << std::endl;
                    merge_toml_tables(v, toml::parse(fname));
                }
            }
            else
            {
                const auto& fname = toml::find<std::string>(v, "include");
                merge_toml_tables(v, toml::parse(fname));
            }
        }
    }
    else if(v.is_array()) // handle an array of tables
    {
        for(auto& elem : v.as_array())
        {
            expand_include(elem);
        }
    }
    return;
}

#endif // OPEN_AICG2_PLUS_INPUT_UTILITY_HPP
