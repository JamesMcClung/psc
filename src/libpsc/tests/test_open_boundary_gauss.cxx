#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc.hxx"
#include "output_fields.hxx"
#include "../psc_config.hxx"
#include "../psc_particle_injectors/dirichlet_boundary_injector.hxx"
#include "../psc_particle_injectors/neumann_boundary_injector.hxx"

// Separate from the injector tests: a binary with more than 4
// integrators crashes at exit (notes/bug-many-integrators-crash.md).

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;

Grid_t* setupGrid()
{
  auto domain = Grid_t::Domain{{1, 8, 2},  // n grid points
                               {1, 8, 2},  // physical lengths
                               {0, 0, 0},  // location of lower corner
                               {1, 1, 1}}; // n patches

  auto bc = psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                          {BND_FLD_PERIODIC, BND_FLD_OPEN, BND_FLD_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC},
                          {BND_PRT_PERIODIC, BND_PRT_OPEN, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(NR_KINDS);
  kinds[KIND_ELECTRON] = {-1.0, 1.0, "e"};
  kinds[KIND_ION] = {1.0, 1.0, "i"};

  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = 1;

  double dt = 1;
  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// ======================================================================
// Gauss's law at the lower boundary node (y = 0)
//
// The gauss check exempts boundary nodes, since injected and deleted particles
// leave virtual charge in the ghost cell. The injector must make that virtual
// charge sit entirely beyond the boundary node, so Gauss's law holds there for
// the real particles alone.

template <typename Mparticles, typename MfieldsState>
double gauss_err_at_lo_node(Mparticles& mprts, MfieldsState& mflds)
{
  const auto& grid = mprts.grid();
  using ItemRho = Checks::item_rho_type;
  auto rho = psc::mflds::interior(grid, ItemRho{grid}(mprts));
  auto dive = psc::mflds::interior(grid, Item_dive<MfieldsState>{}(mflds));
  auto diff = gt::eval(dive - rho);
  return gt::norm_linf(diff.view(_all, _s(0, 1), _all, _all, _all));
}

// Places particles at a given fractional depth into the ghost cell
struct ParticleGeneratorAtDepth
{
  using Real = psc::particle::Inject::Real;
  using Real3 = psc::particle::Inject::Real3;

  psc::particle::Inject get(Real3 min_pos, Real3 pos_range)
  {
    Real uy = n_injected++ < max_n_injected ? 2.0 : 0.0;
    Real3 x = min_pos + pos_range * Real3{.5, 1 - depth, .5};
    return {x, {0.0, uy, 0.0}, 1.0, kind_idx, 0};
  }

  int kind() const { return kind_idx; }

  Real depth;
  int kind_idx;
  int max_n_injected;
  int n_injected = 0;
};

TEST(DirichletBoundaryInjectorTest, GaussAtBoundaryNodeAfterInjection)
{
  PscParams psc_params;
  psc_params.nmax = 1;

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  Checks checks{grid, MPI_COMM_WORLD, ChecksParams{}};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);
  // deep enough that half the charge starts on the ghost node
  psc.add_injector(new DirichletBoundaryInjector<LoHi::Lo, ParticleGeneratorAtDepth,
                                                 PscConfig::PushParticles>(
    ParticleGeneratorAtDepth{.5, KIND_ION, 1}));

  psc.pre_first_step();
  psc.step();

  ASSERT_EQ(mprts.size(), 1);
  EXPECT_LT(gauss_err_at_lo_node(mprts, mflds), 1e-10);
}

TEST(NeumannBoundaryInjectorTest, GaussAtBoundaryNodeAfterCopy)
{
  PscParams psc_params;
  psc_params.nmax = 1;

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};
  Checks checks{grid, MPI_COMM_WORLD, ChecksParams{}};
  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);
  psc.add_injector(
    new NeumannBoundaryInjector<LoHi::Lo, PscConfig::PushParticles>{});

  {
    // a neutral pair; the ion moves from the edge cell inward, so a copy of it
    // is injected that starts in the ghost cell, at y = -.1
    auto injector = mprts.injector();
    auto inj = injector[0];
    inj({{.5, .9, .5}, {0, 0, 0}, 1, KIND_ELECTRON});
    inj({{.5, .9, .5}, {0, 2, 0}, 1, KIND_ION});
  }

  psc.pre_first_step();
  psc.step();

  ASSERT_EQ(mprts.size(), 3);
  EXPECT_LT(gauss_err_at_lo_node(mprts, mflds), 1e-10);
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
