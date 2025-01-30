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

    Topology   topology(7);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_harmonic_bond_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLGaussianBondForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction = "BondLength"
        potential   = "Gaussian"
        topology    = "bond"
        parameters  = [
            {indices = [ 0, 1], k = 1.0, v0 = 2.0, sigma = 7.0},
            {indices = [ 3, 5], k = 3.0, v0 = 4.0, sigma = 8.0},
            {indices = [ 2, 6], k = 5.0, v0 = 6.0, sigma = 9.0}
        ]
    )"_toml;

    Topology   topology(7);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_gaussian_bond_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLGoContactForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction = "BondLength"
        potential   = "GoContact"
        topology    = "bond"
        parameters  = [
            {indices = [ 0, 1], k = 1.0, v0 = 2.0},
            {indices = [ 3, 5], k = 3.0, v0 = 4.0},
            {indices = [ 2, 6], k = 5.0, v0 = 6.0}
        ]
    )"_toml;

    Topology   topology(7);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_go_contact_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLCappedGoContactForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "BondLength"
        potential     = "CappedGoContact"
        topology      = "bond"
        capping_ratio = 0.5
        parameters  = [
            {indices = [ 0, 1], k = 1.0, v0 = 2.0},
            {indices = [ 3, 5], k = 3.0, v0 = 4.0},
            {indices = [ 2, 6], k = 5.0, v0 = 6.0}
        ]
    )"_toml;

    Topology   topology(7);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_capped_go_contact_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLThreeSPN2BondForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "BondLength"
        potential     = "3SPN2Bond"
        topology      = "bond"
        parameters  = [
            {indices = [ 0, 1], k = 1.0, v0 = 2.0},
            {indices = [ 3, 5], k = 3.0, v0 = 4.0},
            {indices = [ 2, 6], k = 5.0, v0 = 6.0}
        ]
    )"_toml;

    Topology   topology(7);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_3spn2_bond_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLHarmonicAngleForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "BondAngle"
        potential     = "Harmonic"
        topology      = "angle"
        parameters  = [
            {indices = [ 0, 1, 2], k = 1.0, v0 = 2.0},
            {indices = [ 3, 4, 5], k = 3.0, v0 = 4.0},
            {indices = [ 2, 6, 7], k = 5.0, v0 = 6.0}
        ]
    )"_toml;

    Topology   topology(8);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_harmonic_angle_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLFlexibleLocalAngleForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "BondAngle"
        potential     = "FlexibleLocalAngle"
        topology      = "none"
        parameters  = [
            {indices = [ 0, 1, 2], k = 1.0, y = "y1_PHE", d2y = "y2_PHE"},
            {indices = [ 3, 4, 5], k = 3.0, y = "y1_VAL", d2y = "y2_VAL"},
            {indices = [ 2, 6, 7], k = 5.0, y = "y1_ALA", d2y = "y2_ALA"}
        ]
    )"_toml;

    Topology topology(8);
    const std::string            aa_type      = "PHE";
    const std::array<double, 10> spline_table{
        5.00, 1.66, 0.92, 0.99, 0.91, 0.97, 1.01, 1.39, 3.23, 10.0};
    const bool                   use_periodic = false;

    EXPECT_NO_THROW(read_toml_flexible_local_angle_ff_generator(
                v, aa_type, spline_table, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLGaussianDihedralForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "DihedralAngle"
        potential     = "Gaussian"
        topology      = "none"
        parameters  = [
            {indices = [ 0,  1,  2,  3], v0 = -1.0, k = -4.0, "σ" = 7.0},
            {indices = [ 1,  2,  3,  4], v0 = -2.0, k = -5.0, "σ" = 8.0},
            {indices = [ 2,  3,  4,  5], v0 = -3.0, k = -6.0, "σ" = 9.0},
        ]
    )"_toml;

    Topology   topology(6);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_gaussian_dihedral_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLCosineDihedralForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "DihedralAngle"
        potential     = "Cosine"
        topology      = "none"
        parameters  = [
            {indices = [ 0,  1,  2,  3], v0 = -1.0, k = -4.0, n = 7},
            {indices = [ 1,  2,  3,  4], v0 = -2.0, k = -5.0, n = 8},
            {indices = [ 2,  3,  4,  5], v0 = -3.0, k = -6.0, n = 9},
        ]
    )"_toml;

    Topology   topology(6);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_cosine_dihedral_ff_generator(v, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLFlexibleLocalDihedralForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction   = "DihedralAngle"
        potential     = "FlexibleLocalDihedral"
        topology      = "none"
        parameters  = [
            {indices = [ 0,  1,  2,  3], k = 1.0, coef = "PHE-VAL"},
            {indices = [ 1,  2,  3,  4], k = 2.0, coef = "VAL-ALA"},
            {indices = [ 2,  3,  4,  5], k = 3.0, coef = "ALA-LEU"},
        ]
    )"_toml;

    Topology                                  topology(6);
    const std::pair<std::string, std::string> aa_pair{"PHE", "R3"};
    const std::array<double, 7 >              fourier_table{
        2.2356,  0.4119, -0.1283,  0.0229, -0.2708, -0.0085, -0.0641};
    const bool                                use_periodic = false;

    EXPECT_NO_THROW(read_toml_flexible_local_dihedral_ff_generator(
                v, aa_pair, fourier_table, topology, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOML3SPN2BaseStackingForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v_3spn2 = u8R"(
        interaction   = "3SPN2BaseStacking"
        potential     = "3SPN2"
        topology      = "nucleotide"
        parameters  = [
            {strand = 0, nucleotide = 0,        S = 0, B = 1, Base = "A"},
            {strand = 0, nucleotide = 1, P = 2, S = 3, B = 4, Base = "T"},
        ]
    )"_toml;

    Topology   topology_3spn2(6);
    const bool use_periodic = false;

    EXPECT_NO_THROW(read_toml_3spn2_base_stacking_ff_generator(
                v_3spn2, topology_3spn2, use_periodic));

    const toml::value v_3spn2c = u8R"(
        interaction   = "3SPN2BaseStacking"
        potential     = "3SPN2C"
        topology      = "nucleotide"
        parameters  = [
            {strand = 0, nucleotide = 0,        S = 0, B = 1, Base = "A"},
            {strand = 0, nucleotide = 1, P = 2, S = 3, B = 4, Base = "T"},
        ]
    )"_toml;

    Topology topology_3spn2c(6);

    EXPECT_NO_THROW(read_toml_3spn2_base_stacking_ff_generator(
                v_3spn2c, topology_3spn2c, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLExcludedVolumeForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction                     = "Pair"
        potential                       = "ExcludedVolume"
        ignore.molecule                 = "Nothing"
        ignore.particles_within.bond    = 3
        ignore.particles_within.contact = 1
        epsilon                         = 0.6
        cutoff                          = 0.5
        parameters = [
        {index = 0, radius = 1.0},
        {index = 1, radius = 2.0},
        {index = 2, radius = 3.0},
        {index = 3, radius = 4.0},
        {index = 4, radius = 5.0},
        {index = 5, radius = 6.0},
        {index = 6, radius = 7.0},
        ]
    )"_toml;

    const std::size_t system_size(7);
    Topology          topology(system_size);
    topology.add_edges({{0, 1}, {1, 2}, {2, 3}, {3, 4}, {5, 6}}, "bond");
    const std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    const bool                                    use_periodic = false;

    EXPECT_NO_THROW(read_toml_excluded_volume_ff_generator(
                v, system_size, topology, group_vec, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOML3SPN2ExcludedVolumeForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction     = "Pair"
        potential       = "3SPN2ExcludedVolume"
        ignore.molecule = "Nothing"
        ignore.particle_within.bond       = 1
        ignore.particle_within.nucleotide = 1
        parameters = [
        # strand 1
        {index = 0, kind = "S"},
        {index = 1, kind = "A"},
        {index = 2, kind = "P"},
        # strand 2
        {index = 3, kind = "S"},
        {index = 4, kind = "T"},
        {index = 5, kind = "P"},
        {index = 6, kind = "S"},
        {index = 7, kind = "A"},
        {index = 8, kind = "P"},
        ]
    )"_toml;

    const std::size_t system_size(9);
    Topology          topology(system_size);
    topology.add_edges({{0, 1}, {1, 2},
                        {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}}, "bond");
    std::vector<std::pair<std::size_t, std::size_t>> edges;
    edges.push_back(std::make_pair(0, 3));
    edges.push_back(std::make_pair(0, 4));
    edges.push_back(std::make_pair(1, 3));
    edges.push_back(std::make_pair(1, 4));
    edges.push_back(std::make_pair(2, 3));
    edges.push_back(std::make_pair(2, 4));
    edges.push_back(std::make_pair(0, 5));
    edges.push_back(std::make_pair(1, 5));
    edges.push_back(std::make_pair(2, 5));
    topology.add_edges(edges, "nucleotide");

    const std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    const bool                                    use_periodic = false;

    EXPECT_NO_THROW(read_toml_3spn2_excluded_volume_ff_generator(
                v, system_size, topology, group_vec, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLWeeksChandlerAndersenForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction             = "Pair"
        potential               = "WCA"
        ignore_particles_within = {bond = 1, angle = 1}
        env.test_epsilon = 0.5
        env.test_sigma_A = 7.0
        env.test_sigma_B = 4.0
        parameters = [
            {index = 0, sigma = "test_sigma_A", epsilon = "test_epsilon"},
            {index = 1, sigma = "test_sigma_B", epsilon = "test_epsilon"},
            {index = 2, sigma = "test_sigma_A", epsilon = "test_epsilon"},
            {index = 3, sigma = "test_sigma_B", epsilon = "test_epsilon"}
        ]
    )"_toml;

    const std::size_t system_size(4);
    Topology          topology(system_size);
    topology.add_edges({{0, 1}, {1, 2}, {2, 3}}, "bond");
    topology.add_edges({{0, 1}, {1, 2}, {2, 3}}, "angle");

    const std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    const bool                                    use_periodic = false;

    EXPECT_NO_THROW(read_toml_weeks_chandler_andersen_ff_generator(
                v, system_size, topology, group_vec, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLUniformWeeksChandlerAndersenForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction = "Pair"
        potential   = "WCA"
        table.TESTA.TESTB = {sigma = 5.0, epsilon = 0.5}
        parameters = [
            {index = 0, name = "TESTA"},
            {index = 1, name = "TESTB"},
            {index = 2, name = "TESTA"},
            {index = 3, name = "TESTB"}
        ]
    )"_toml;

    const std::size_t system_size(4);
    const bool        use_periodic = false;
    Topology          topology(system_size);
    const std::pair<std::string, std::string>&    name_pair{"TESTA", "TESTB"};
    const std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);
    const std::vector<std::pair<std::size_t, std::size_t>> ignore_list{};
    const std::vector<std::pair<std::string, std::string>> ignore_group_pairs{};

    EXPECT_NO_THROW(read_toml_uniform_weeks_chandler_andersen_ff_generator(
                v, system_size, 5.0/*sigma*/, 0.5/*epsilon*/,
                name_pair, topology, ignore_list, ignore_group_pairs,
                group_vec, use_periodic));
}

TEST(ReadTOMLForceFieldGenerator, ReadTOMLDebyeHuckelForceFieldGenerator)
{
    using namespace toml::literals;
    const toml::value v = u8R"(
        interaction = "Pair"
        potential   = "WCA"
        ignore.particles_within.bond = 3
        parameters = [
        {index = 2, charge = -1.0},
        {index = 5, charge = -1.0},
        {index = 8, charge = -1.0}
        ]
    )"_toml;

    const std::size_t system_size(9);
    const bool        use_periodic = false;
    Topology          topology(system_size);
    topology.add_edges({{0, 1}, {1, 2}, {3, 4}, {4, 5},
                        {5, 6}, {6, 7}, {7, 8}}, "bond");
    const double      ionic_strength(0.15);
    const double      temperature(300.0);
    const std::vector<std::optional<std::string>> group_vec(system_size, std::nullopt);

    EXPECT_NO_THROW(read_toml_debye_huckel_ff_generator(
                v, system_size, ionic_strength, temperature, topology,
                group_vec, use_periodic));
}
