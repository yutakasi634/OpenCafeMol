#ifndef OPEN_AICG2_PLUS_PROGRESS_BAR_HPP
#define OPEN_AICG2_PLUS_PROGRESS_BAR_HPP

#include <array>
#include <chrono>
#include <ostream>

#include <cmath>
#include <cstdio>

class ProgressBar
{
  public:

    ProgressBar(std::size_t bar_width) noexcept
        : bar_width_(bar_width), start_(std::chrono::system_clock::now())
    {}

    void format(std::size_t count, std::size_t total_step, std::ostream& os) const
    {
        const double ratio = (count == total_step) ? 1.0 :
            std::max(0.0, std::min(1.0, count / static_cast<double>(total_step)));

        std::array<char, 8> buf;
        buf.fill('\0');
        std::snprintf(buf.data(), 8, "%5.1f%%|", ratio * 100.0);
        os << '\r' << buf.data();

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
            os << full;
        }
        if(filled < bar_width_)
        {
            switch(static_cast<std::size_t>(ratio * bar_width_ * 8) % 8)
            {
                case 0: {os << ' '; break;}
                case 1: {os << l1;  break;}
                case 2: {os << l2;  break;}
                case 3: {os << l3;  break;}
                case 4: {os << l4;  break;}
                case 5: {os << l5;  break;}
                case 6: {os << l6;  break;}
                case 7: {os << l7;  break;}
            }
            for(std::size_t i=1; i<bar_width_-filled; ++i)
            {
                os << ' ';
            }
        }
        os << '|';

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
                os << buf.data() << " seconds remaining  ";
            }
            else if(residual < 60.0 * 60.0)
            {
                std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167);
                os << buf.data() << " minutes remaining ";
            }
            else if(residual < 60.0 * 60.0 * 24.0)
            {
                std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167 * 0.0167);
                os << buf.data() << " hours remaining   ";
            }
            else if(residual < 60.0 * 60.0 * 24.0 * 99.0)
            {
                std::snprintf(buf.data(), 6, "%5.1f", residual * 0.0167 * 0.0167 * 0.0417);
                os << buf.data() << " days remaining    ";
            }
            else
            {
                os << " over 100 days remaining";
            }
        }
        os << std::flush;

        return;
    }

    void finalize(std::ostream& os) const
    {
        this->format(1, 1, os); // show 100%
        os << std::endl;
    }

  private:
    std::size_t                           bar_width_;
    std::chrono::system_clock::time_point start_;
};

#endif // OPEN_AICG2_PLUS_PROGRESS_BAR_HPP
