#ifndef OPEN_AICG2_PLUS_FORCE_FIELD_ID_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FORCE_FIELD_ID_GENERATOR_HPP

#include <atomic>
#include <cstddef>

class ForceFieldIDGenerator
{
  public:

    ForceFieldIDGenerator() noexcept
        : id_(0)
    {}

    std::size_t gen() noexcept
    {
        return id_++;
    }

  private:

    // to make it thread-safe (maybe we don't need to make it thread safe, though)
    std::atomic<std::size_t> id_;
};

inline ForceFieldIDGenerator ffid;

#endif// OPEN_AICG2_PLUS_FORCE_FIELD_ID_GENERATOR_HPP
