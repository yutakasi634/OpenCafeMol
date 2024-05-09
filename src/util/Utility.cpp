#include "Utility.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>

namespace Utility
{

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
    std::filesystem::path fpath(filename);
    fpath.remove_filename(); // extract direcotry

    // output filename = "filename.pdb" means output dir is current dir
    if(fpath.empty())
    {
        fpath = "./";
    }

    if( ! std::filesystem::exists(fpath))
    {
        std::filesystem::create_directories(fpath);
    }
    if( ! std::filesystem::exists(fpath) || ! std::filesystem::is_directory(fpath))
    {
        throw std::runtime_error(
                "failed to make output directory: " + fpath.string());
    }

    std::ofstream ofs(filename);
    if(not ofs.good())
    {
        throw std::runtime_error("file open error : " + filename);
    }
    ofs.close();
    return;
}

const toml::value& find_either(const toml::value& v,
        const std::string& key1, const std::string& key2)
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


} // Utility
