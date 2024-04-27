#ifndef OPEN_AICG2_PLUS_PROGRESS_BAR_HPP
#define OPEN_AICG2_PLUS_PROGRESS_BAR_HPP

#include <chrono>


class ProgressBar
{
  public:

    ProgressBar(std::size_t bar_width) noexcept
        : bar_width_(bar_width), start_(std::chrono::system_clock::now())
    {}

    void format(std::size_t count, std::size_t total_step) const;

    void finalize() const;

  private:
    std::size_t                           bar_width_;
    std::chrono::system_clock::time_point start_;
};

#endif // OPEN_AICG2_PLUS_PROGRESS_BAR_HPP
