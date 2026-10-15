#pragma once

#include <gtest/gtest.h>

#include "test_common.hxx"

#include "psc.hxx"
#include "../psc_config.hxx"
#include "../psc_particle_injectors/neumann_boundary_injector.hxx"

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;

Grid_t* setupGrid(double cfl)
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

  double dt = cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// For each (y, uy), places a stationary electron and a moving ion at the same
// position (to satisfy Gauss' law at t=0), runs the simulation for `nmax` steps
// with a NeumannBoundaryInjector at the given boundary, and returns the
// final particles' y positions, sorted. Copies are made at `copy_depth`.
template <LoHi LOHI, typename Resampler = NeumannResamplerNone>
std::vector<double> run(std::vector<std::pair<double, double>> ys_uys, int nmax,
                        Resampler resampler = {}, int copy_depth = 1)
{
  PscParams psc_params;
  psc_params.nmax = nmax;
  psc_params.stats_every = 1;
  psc_params.cfl = .75;

  auto grid_ptr = setupGrid(psc_params.cfl);
  auto& grid = *grid_ptr;

  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  ChecksParams checks_params{};
  checks_params.continuity.check_interval = 1;
  checks_params.gauss.check_interval = 1;
  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  Balance balance{.1};
  Collision collision{grid, 0, 0.1};

  auto psc = makePscIntegrator<PscConfig>(psc_params, grid, mflds, mprts,
                                          balance, collision, checks);

  NeumannBoundaryInjector<LOHI, PscConfig::PushParticles, Resampler> injector{
    resampler, copy_depth};
  psc.add_injector(&injector);

  {
    auto inj = mprts.injector()[0];
    for (auto [y, uy] : ys_uys) {
      inj({{0, y, .5}, {0, 0, 0}, 1, KIND_ELECTRON});
      inj({{0, y, .5}, {0, uy, 0}, 1, KIND_ION});
    }
  }

  psc.pre_first_step();
  for (; grid.timestep_ < psc_params.nmax;) {
    psc.step();

    EXPECT_LT(checks.continuity.last_max_err, checks.continuity.err_threshold);
    EXPECT_LT(checks.gauss.last_max_err, checks.gauss.err_threshold);
  }

  std::vector<double> ys;
  auto accessor = mprts.accessor();
  for (auto prt : accessor[0]) {
    ys.push_back(prt.position()[1]);
  }
  std::sort(ys.begin(), ys.end());
  return ys;
}
