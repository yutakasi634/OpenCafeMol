#include "ProgressBar.hpp"

#include <array>
#include <iostream>

#include <cstdio>
#include <cmath>

void ProgressBar::format(std::size_t count, std::size_t total_step) const
{
    const double ratio = (count == total_step) ? 1.0 :
        std::max(0.0, std::min(1.0, count / static_cast<double>(total_step)));

    std::array<char, 8> buf;
    buf.fill('\0');
    std::snprintf(buf.data(), 8, "%5.1f%%|", ratio * 100.0);
    std::cerr << '\r' << buf.data();

    constexpr auto full  = u8"█"; // U+2588 Full block
    constexpr auto l7    = u8"▉"; // U+2589 Left seven eighths block
    constexpr auto l6    = u8"▊"; // U+258A Left three quarters block
    constexpr auto l5    = u8"▋"; // U+258B Left five eighths block
    constexpr auto l4    = u8"▌"; // U+258C Left half block
    constexpr auto l3    = u8"▍"; // U+258D Left three eighths block
    constexpr auto l2    = u8"▎"; // U+258E Left one quarter block
    constexpr auto l1    = u8"▏"; // U+258F Left one eighth block

    const std::size_t filled = std::floor(ratio * bar_width_);
    for(std::size_t i=0; i<filled; ++i)
    {
        std::cerr << full;
    }
    if(filled < bar_width_)
    {
        switch(static_cast<std::size_t>(ratio * bar_width_ * 8) % 8)
        {
            case 0: {std::cerr << ' '; break;}
            case 1: {std::cerr << l1;  break;}
            case 2: {std::cerr << l2;  break;}
            case 3: {std::cerr << l3;  break;}
            case 4: {std::cerr << l4;  break;}
            case 5: {std::cerr << l5;  break;}
            case 6: {std::cerr << l6;  break;}
            case 7: {std::cerr << l7;  break;}
        }
        for(std::size_t i=1; i<bar_width_-filled; ++i)
        {
            std::cerr << ' ';
        }
    }
    std::cerr << '|';

    const auto current  = std::chrono::system_clock::now();

    buf.fill('\0');

    // ratio == 0 means the first step.
    // If ratio is zero, we cannot estimate how long does it take.
    // Also we need to avoid zero-division.
    if(ratio != 0)
    {
        const auto residual = std::chrono::duration_cast<std::chrono::milliseconds>(
            (current - start_) * (1.0 - ratio) / ratio).count() * 0.001;
        if(residual < 60.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual);
            std::cerr << buf.data() << " seconds remaining  ";
        }
        else if(residual < 60.0 * 60.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167);
            std::cerr << buf.data() << " minutes remaining ";
        }
        else if(residual < 60.0 * 60.0 * 24.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167 * 0.0167);
            std::cerr << buf.data() << " hours remaining   ";
        }
        else if(residual < 60.0 * 60.0 * 24.0 * 99.0)
        {
            std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167 * 0.0167 * 0.0417);
            std::cerr << buf.data() << " days remaining    ";
        }
        else
        {
            std::cerr << " over 100 days remaining";
        }
    }
    std::cerr << std::flush;

    return;
}

void ProgressBar::finalize() const
{
    this->format(1, 1); // show 100%
    std::cerr << std::endl;
}

