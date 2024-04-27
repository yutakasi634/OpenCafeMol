#include "Logger.hpp"

#include <iostream>

namespace detail
{

void log_impl(std::string log)
{
    std::cerr << log << std::endl;
    return ;
}

} // detail
