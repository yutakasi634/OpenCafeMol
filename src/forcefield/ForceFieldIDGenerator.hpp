#ifndef OPEN_AICG2_PLUS_FORCE_FIELD_ID_GENERATOR_HPP
#define OPEN_AICG2_PLUS_FORCE_FIELD_ID_GENERATOR_HPP

#include <atomic>
#include <cstddef>

// OpenMM does not introduce namespace to parameters defined in CustomForces.
// So we need to use Hungarian notation to parameters to avoid collisions.
// Especially, since we define multiple FLA/FLD Forces, we need to add
// forcefield-unique index to this notation, not only the forcefield name.
// On some platform, OpenMM does not check parameter collision (CUDA does not,
// but CPU does) before running simulation. It causes a fatal but silent bug.
// To avoid this, we use ID generator that automatically generates unique idx.
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
