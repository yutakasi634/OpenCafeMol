#ifndef OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP
#define OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP

#include "src/util/Constants.hpp"
#include "ForceFieldGeneratorBase.hpp"
#include "ForceFieldIDGenerator.hpp"

#include <OpenMM.h>
#include <fmt/core.h>

#include <iostream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>


class ProteinDNANonSpecificForceFieldGenerator final : public ForceFieldGeneratorBase
{
  public:
    using indices_dna_type     = std::array<std::size_t, 2>;
    using indices_protein_type = std::array<std::size_t, 3>;

  public:
    ProteinDNANonSpecificForceFieldGenerator(
        const std::vector<indices_dna_type>& indices_dna,
        const std::vector<indices_protein_type>& indices_protein,
        const double sigma, const double delta, const double cutoff_ratio,
        const std::vector<double> ks, const std::vector<double> r0s,
        const std::vector<double> theta0s, const std::vector<double> phi0s,
        const bool use_periodic);

    std::unique_ptr<OpenMM::Force> generate() const noexcept override;

    std::string name() const noexcept { return "PDNS"; }

  private:
    std::vector<indices_dna_type> indices_dna_;
    std::vector<indices_protein_type> indices_protein_;
    double sigma_;
    double delta_;
    double cutoff_ratio_;
    std::vector<double> ks_;
    std::vector<double> r0s_;
    std::vector<double> theta0s_;
    std::vector<double> phi0s_;
    bool                use_periodic_;
    std::string         ffgen_id_;
};

#endif // OPEN_AICG2_PLUS_PROTEIN_DNA_NON_SPECIFIC_FORCE_FIELD_GENERATOR_HPP
