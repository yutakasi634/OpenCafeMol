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

inline void
add_offset(std::size_t& index, const toml::value& offset)
{
    if(offset.is_array())
    {
        throw std::runtime_error("[error] add_offset: "
            "offset array is invalid here. expected type is scalar.");
    }
    else
    {
        index += offset.as_integer();
    }
}

inline void
add_offset(std::pair<std::size_t, std::size_t>& indices, const toml::value& offset)
{
    if(offset.is_array())
    {
        if(offset.size() != 2)
        {
            throw std::runtime_error("[error] add_offset: "
                "invalid size of offset array " + std::to_string(offset.size()) +
                 " here. the size of offset array must much to the size of "
                 "indices array. expected size is 2.");
        }

        indices.first  += offset[0].as_integer();
        indices.second += offset[1].as_integer();
    }
    else
    {
        const std::size_t offset_val = offset.as_integer();
        indices.first  += offset_val;
        indices.second += offset_val;
    }
}

template <std::size_t N>
void add_offset(std::array<std::size_t, N>& indices, const toml::value& offset)
{
    if(offset.is_array())
    {
        if(offset.size() != N)
        {
            throw std::runtime_error("[error] add_offset: "
                "invalid size of offset array " + std::to_string(offset.size()) +
                 " here. the size of offset array must much to the size of "
                 "indices array. expected size is " +
                 std::to_string(indices.size()) + ".");
        }

        for(std::size_t i=0; i<N; ++i)
        {
            indices[i] += offset[i].as_integer();
        }
    }
    else
    {
        for(auto& i : indices) {i += offset.as_integer();}
    }
}

#endif // OPEN_AICG2_PLUS_INPUT_UTILITY_HPP
