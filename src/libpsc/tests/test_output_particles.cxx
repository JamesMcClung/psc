
#include "gtest/gtest.h"

#include <dim.hxx>
#include <psc_particles_single.h>
#include <output_particles.hxx>

struct OutpAscii
{
  static constexpr const char* type = "ascii";
};

struct OutpHdf5Single
{
  static constexpr const char* type = "hdf5_single";
};

struct OutpHdf5Double
{
  static constexpr const char* type = "hdf5_double";
};

template<typename _Dim, typename _Mparticles, typename _Outp>
struct Config
{
  using Dim = _Dim;
  using Mparticles = _Mparticles;
  using Outp = _Outp;
};

template<typename T>
struct OutputParticlesTest : ::testing::Test
{
  using Dim = typename T::Dim;

  const double L = 160;

  void make_psc(const Grid_t::Kinds& kinds)
  {
    Int3 gdims = {16, 16, 16};
    if (Dim::InvarX::value) { gdims[0] = 1; ibn[0] = 0; }
    if (Dim::InvarY::value) { gdims[1] = 1; ibn[1] = 0; }
    if (Dim::InvarZ::value) { gdims[2] = 1; ibn[2] = 0; }

    auto grid_domain = Grid_t::Domain{gdims, {L, L, L}};
    auto grid_bc = GridBc{{ BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC },
			  { BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC }};
    
    auto norm_params = Grid_t::NormalizationParams::dimensionless();
    norm_params.nicell = 200;
    auto coeff = Grid_t::Normalization{norm_params};

    grid_ = Grid_t::psc_make_grid(grid_domain, grid_bc, kinds, coeff, 1., {});
    ggrid = grid_; // FIXME
  }

  const Grid_t& grid() const
  {
    assert(grid_);
    return *grid_;
  }

private:
  Grid_t* grid_;
  Int3 ibn = { 2, 2, 2 };
};

using OutputParticlesTestTypes = ::testing::Types<Config<dim_xyz, MparticlesSingle, OutpAscii>,
						  Config<dim_xyz, MparticlesSingle, OutpHdf5Single>,
						  Config<dim_xyz, MparticlesSingle, OutpHdf5Double>>;

TYPED_TEST_CASE(OutputParticlesTest, OutputParticlesTestTypes);

// ======================================================================
// Test1

TYPED_TEST(OutputParticlesTest, Test1)
{
  using Mparticles = typename TypeParam::Mparticles;

  auto kinds = Grid_t::Kinds{{1., 100., "ion"}, {-1., 1., "electron"}};
  this->make_psc(kinds);
  const auto& grid = this->grid();
  
  // init particle
  auto n_prts_by_patch = std::vector<uint>{1};

  Mparticles mprts{grid};
  {
    auto injector = mprts[0].injector();
    injector({{1., 0., 0.}, {}, 1., 0});
    injector({{2., 0., 0.}, {}, 1., 1});
  }

  psc_output_particles* outp = psc_output_particles_create(grid.comm());
  psc_output_particles_set_type(outp, TypeParam::Outp::type);
  psc_output_particles_set_param_int(outp, "every_step", 1);
  psc_output_particles_setup(outp);
  //psc_output_particles_view(outp);
  PscOutputParticlesBase{outp}.run(mprts);
  psc_output_particles_destroy(outp);
}

// ----------------------------------------------------------------------
// main

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}


