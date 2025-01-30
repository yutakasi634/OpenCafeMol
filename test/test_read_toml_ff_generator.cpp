#include "src/input/ReadTOMLForceFieldGenerator.hpp"

#include <toml.hpp>

#include <gtest/gtest.h>


TEST(ReadTOMLForceFieldGenerator, ReadTOMLHarmonicBondForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction = "BondLength"
        potential   = "Harmonic"
        topology    = "bond"
        parameters  = [
            {indices = [ 0, 1], v0 = 1.0, k = 2.0},
            {indices = [ 3, 5], v0 = 3.0, k = 4.0},
            {indices = [ 2, 6], v0 = 5.0, k = 6.0}
        ]
    )"_toml;

    Topology topology(7);
    bool     use_periodic = false;

    EXPECT_NO_THROW(read_toml_harmonic_bond_ff_generator(v, topology, use_periodic));
}
