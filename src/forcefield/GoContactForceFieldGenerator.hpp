#ifndef OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP

#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

class GoContactForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_type = std::pair<std::size_t, std::size_t>;

  public:
    GoContactForceFieldGenerator(
        const std::vector<indices_type>& indices_vec,
        const std::vector<double>& ks, const std::vector<double>& r0s,
        const bool use_periodic)
        : indices_vec_(indices_vec), ks_(ks), r0s_(r0s),
          use_periodic_(use_periodic), ffgen_id_(fmt::format("GC{}", ffid.gen()))
    {
        if(!(indices_vec.size() == ks.size() && ks.size() == r0s.size()))
        {
            std::ostringstream oss;
            oss << "[error] GoContactForceFieldGenerator : "
                   "parameter number of "
                   "indices_vec (" << indices_vec.size() << "), "
                   "k ("           << ks.size()          << ") and "
                   "r0 ("          << r0s.size()         << ") is not matched."
                   "The number of these parameters must be same.";
            throw std::runtime_error(oss.str());
        }
    }

    std::unique_ptr<OpenMM::Force> generate() const override;

    const std::vector<indices_type>& indices() const noexcept { return indices_vec_; }
    std::string name() const override { return "GoContact"; }

  private:
    std::vector<indices_type> indices_vec_;
    std::vector<double>       ks_;
    std::vector<double>       r0s_;
    bool                      use_periodic_;
    std::string               ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_GOCONTACT_FORCE_FIELD_GENERATOR_HPP
