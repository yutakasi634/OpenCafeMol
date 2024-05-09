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

    void format(std::size_t count, std::size_t total_step, std::ostream& os) const;

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
