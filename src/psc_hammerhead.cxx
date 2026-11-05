#include <psc.hxx>
#include <setup_fields.hxx>
#include <setup_particles.hxx>

#include "output_fields.hxx"
#include "psc_config.hxx"
#include "input_params.hxx"
#include "libpsc/psc_output_particles/output_particles_adios2_impl.hxx"

// ======================================================================
// psc_hammerhead
//
// Periodic 2d box in the y-z plane with a uniform background magnetic field
// along y. The plasma consists of thermal electrons and an ion core + ion
// beam, both bi-Maxwellian with respect to B. The beam drifts along B; the
// core drifts the opposite way such that the net ion drift is 0. Electrons
// don't drift, so there is no net current.
//
// Densities: n_e = n_core + n_beam = 1.

// ======================================================================
// PSC configuration
//
// This sets up compile-time configuration for the code, in particular
// what data structures and algorithms to use
//
// EDIT to change order / floating point type / cuda / 2d/3d

using Dim = dim_yz;
#ifdef USE_CUDA
using PscConfig = PscConfig1vbecCuda<Dim>;
#else
using PscConfig = PscConfig1vbecDouble<Dim>;
#endif

// ----------------------------------------------------------------------

using MfieldsState = PscConfig::MfieldsState;
using Mparticles = PscConfig::Mparticles;
using Balance = PscConfig::Balance;
using Collision = PscConfig::Collision;
using Checks = PscConfig::Checks;
using Marder = PscConfig::Marder;
using OutputParticles = PscConfig::OutputParticles;
using real_t = PscConfig::Mfields::real_t;
using Real3 = Vec3<real_t>;

enum
{
  MY_ELECTRON,
  MY_ION_CORE,
  MY_ION_BEAM,
  N_MY_KINDS,
};

// ======================================================================
// Global parameters

PscParams psc_params;

double electron_mass;
double ion_mass;

double b0; // background field, along y

double te; // isotropic electron temperature

double n_core;
double v_core; // parallel to B, fixed by zero net ion drift
double t_core_par;
double t_core_perp;

double n_beam;
double v_beam; // parallel to B, in the frame with zero net ion drift
double t_beam_par;
double t_beam_perp;

Int3 gdims;
Double3 lengths;
Int3 n_patches;

int out_interval;
int marder_interval;

int nicell;

std::string checkpoint_filename;

// ======================================================================
// setupParameters

void setupParameters(int argc, char** argv)
{
  if (argc != 2) {
    LOG_ERROR("Usage: %s path/to/params.txt\n", argv[0]);
  }
  std::string path_to_params(argv[1]);
  InputParams inputParams(path_to_params);

  psc_params.stats_every =
    inputParams.getOrDefault<int>("stats_interval", 1000);
  psc_params.cfl = inputParams.getOrDefault<double>("cfl", .75);
  checkpoint_filename =
    inputParams.getOrDefault<std::string>("checkpoint_filename", "");

  electron_mass = inputParams.get<double>("m_e");
  ion_mass = inputParams.get<double>("m_i");

  b0 = inputParams.get<double>("B_0");

  te = inputParams.get<double>("T_e");

  n_beam = inputParams.get<double>("n_beam");
  if (n_beam < 0.0 || n_beam >= 1.0) {
    LOG_ERROR("n_beam must be in [0, 1); got %f\n", n_beam);
  }
  v_beam = inputParams.get<double>("v_beam");
  t_beam_par = inputParams.get<double>("T_beam_par");
  t_beam_perp = inputParams.get<double>("T_beam_perp");

  n_core = 1.0 - n_beam;
  v_core = -n_beam * v_beam / n_core;
  t_core_par = inputParams.get<double>("T_core_par");
  t_core_perp = inputParams.get<double>("T_core_perp");

  gdims[0] = inputParams.get<int>("nx");
  gdims[1] = inputParams.get<int>("ny");
  gdims[2] = inputParams.get<int>("nz");
  if (gdims[0] != 1) {
    LOG_ERROR("nx must be 1 (the box is in the y-z plane); got %d\n", gdims[0]);
  }
  psc_params.nmax = inputParams.get<int>("nt");

  n_patches[0] = inputParams.get<int>("npx");
  n_patches[1] = inputParams.get<int>("npy");
  n_patches[2] = inputParams.get<int>("npz");

  if (inputParams.has("lx")) {
    lengths[0] = inputParams.get<double>("lx");
  } else {
    lengths[0] = inputParams.get<double>("dx") * gdims[0];
  }
  if (inputParams.has("ly")) {
    lengths[1] = inputParams.get<double>("ly");
  } else {
    lengths[1] = inputParams.get<double>("dy") * gdims[1];
  }
  if (inputParams.has("lz")) {
    lengths[2] = inputParams.get<double>("lz");
  } else {
    lengths[2] = inputParams.get<double>("dz") * gdims[2];
  }

  if (inputParams.has("checkpoint_interval")) {
    psc_params.write_checkpoint_every_step =
      inputParams.get<int>("checkpoint_interval");
    inputParams.errIfPresent(
      "n_checkpoints",
      "n_checkpoints is mutually exclusive with checkpoint_interval");
  } else if (inputParams.has("n_checkpoints")) {
    int n_checkpoints = inputParams.get<int>("n_checkpoints");
    if (n_checkpoints > 0) {
      psc_params.write_checkpoint_every_step = psc_params.nmax / n_checkpoints;
    }
  }

  int n_writes = inputParams.getOrDefault<int>("n_writes", 100);
  out_interval = psc_params.nmax / n_writes;
  marder_interval = inputParams.getOrDefault<int>("marder_interval", -1);

  nicell = inputParams.getOrDefault<int>("nicell", 100);

  std::ifstream src(path_to_params, std::ios::binary);
  std::ofstream dst("params_record.txt", std::ios::binary);
  dst << src.rdbuf();
}

// ======================================================================
// setupGrid
//
// This helper function is responsible for setting up the "Grid",
// which is really more than just the domain and its decomposition, it
// also encompasses PC normalization parameters, information about the
// particle kinds, etc.

Grid_t* setupGrid()
{
  Double3 corner = {0.0, 0.0, 0.0};
  auto domain = Grid_t::Domain{gdims, lengths, corner, n_patches};

  auto bc =
    psc::grid::BC{{BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_FLD_PERIODIC, BND_FLD_PERIODIC, BND_FLD_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC},
                  {BND_PRT_PERIODIC, BND_PRT_PERIODIC, BND_PRT_PERIODIC}};

  auto kinds = Grid_t::Kinds(N_MY_KINDS);
  kinds[MY_ELECTRON] = {-1.0, electron_mass, "e"};
  kinds[MY_ION_CORE] = {1.0, ion_mass, "core"};
  kinds[MY_ION_BEAM] = {1.0, ion_mass, "beam"};

  // --- generic setup
  auto norm_params = Grid_t::NormalizationParams::dimensionless();
  norm_params.nicell = nicell;

  double dt = psc_params.cfl * courant_length(domain);
  Grid_t::Normalization norm{norm_params};

  int n_ghosts = 2;
  Int3 ibn = n_ghosts * Dim::get_noninvariant_mask();

  return new Grid_t{domain, bc, kinds, norm, dt, -1, ibn};
}

// ======================================================================
// initializeParticles

void initializeParticles(Balance& balance, Grid_t*& grid_ptr, Mparticles& mprts)
{
  SetupParticles<Mparticles> setup_particles(*grid_ptr);
  setup_particles.centerer = centering::Centerer(centering::CC);
  setup_particles.initial_momentum_gamma_correction = true;
  // Keep random_offsets off: every particle then starts at its cell center, so
  // the initial charge density is exactly 0. With random offsets, each kind
  // draws its positions from an identically-seeded stream, which only cancels
  // when all kinds have the same number of particles per cell.

  auto init_np = [&](int kind, Double3 pos, int p, Int3 idx,
                     psc_particle_np& np) {
    // B is along y, so T = {T_perp, T_par, T_perp}
    switch (kind) {
      case MY_ELECTRON:
        np.n = 1.0;
        np.p = setup_particles.createMaxwellian(
          {np.kind, np.n, {0.0, 0.0, 0.0}, {te, te, te}, np.tag});
        break;
      case MY_ION_CORE:
        np.n = n_core;
        np.p = setup_particles.createMaxwellian(
          {np.kind,
           np.n,
           {0.0, v_core, 0.0},
           {t_core_perp, t_core_par, t_core_perp},
           np.tag});
        break;
      case MY_ION_BEAM:
        np.n = n_beam;
        np.p = setup_particles.createMaxwellian(
          {np.kind,
           np.n,
           {0.0, v_beam, 0.0},
           {t_beam_perp, t_beam_par, t_beam_perp},
           np.tag});
        break;
      default: assert(0);
    }
  };

  partitionAndSetupParticles(setup_particles, balance, grid_ptr, mprts,
                             init_np);
}

// ======================================================================
// initializeFields

void add_background_fields(MfieldsState& mflds)
{
  const auto& grid = mflds.grid();

  for (int p = 0; p < mflds.n_patches(); ++p) {
    auto field_patch = make_Fields3d<dim_xyz>(mflds[p]);

    int n_ghosts = mflds.ibn().max();
    grid.Foreach_3d(n_ghosts, n_ghosts, [&](int jx, int jy, int jz) {
      field_patch(HY, jx, jy, jz) += b0;
    });
  }
}

// ======================================================================
// run

static void run(int argc, char** argv)
{
  mpi_printf(MPI_COMM_WORLD, "*** Setting up...\n");

  // ----------------------------------------------------------------------
  // setup various parameters first

  setupParameters(argc, argv);

  // ----------------------------------------------------------------------
  // Set up grid, state fields, particles

  auto grid_ptr = setupGrid();
  auto& grid = *grid_ptr;
  MfieldsState mflds{grid};
  Mparticles mprts{grid};

  // ----------------------------------------------------------------------
  // Set up various objects needed to run this case

  // -- Balance
  psc_params.balance_interval = 1000;
  Balance balance{.1};

  // -- Sort
  psc_params.sort_interval = 100;

  // -- Collision
  int collision_interval = 0;
  double collision_nu = .1;
  Collision collision{grid, collision_interval, collision_nu};

  // -- Checks
  ChecksParams checks_params{};
  checks_params.gauss.check_interval = out_interval;
  checks_params.continuity.check_interval = out_interval;

  Checks checks{grid, MPI_COMM_WORLD, checks_params};

  // -- Marder correction
  double marder_diffusion = 0.9;
  int marder_loop = 3;
  bool marder_dump = false;
  psc_params.marder_interval = marder_interval;
  Marder marder(grid, marder_diffusion, marder_loop, marder_dump);

  // ----------------------------------------------------------------------
  // Set up output

  // -- output fields
  OutputFields<MfieldsState, Mparticles> out_fields;
  out_fields.pfield.out_interval = out_interval;

  OutputMoments<MfieldsState, Mparticles, Dim> out_moments;
  out_moments.pfield.out_interval = out_interval;

  // -- output particles
  OutputParticlesAdios2Params outp_params{};
  outp_params.every_step = out_interval;
  OutputParticlesAdios2<Mparticles, float> outp{grid, outp_params};

  int oute_interval = -100;
  DiagEnergies<Mparticles, MfieldsState> oute{grid.comm(), oute_interval};

  // ----------------------------------------------------------------------
  // set up initial conditions

  if (checkpoint_filename.empty()) {
    initializeParticles(balance, grid_ptr, mprts);
    add_background_fields(mflds);
  } else {
    read_checkpoint(checkpoint_filename, *grid_ptr, mprts, mflds);
  }

  // ----------------------------------------------------------------------
  // run the simulation

  auto psc = makePscIntegrator<PscConfig>(psc_params, *grid_ptr, mflds, mprts,
                                          balance, collision, checks);

  psc.add_gauss_corrector(&marder);

  psc.add_diagnostic(&out_fields);
  psc.add_diagnostic(&out_moments);
  psc.add_diagnostic(&outp);
  psc.add_diagnostic(&oute);

  psc.integrate();
}

// ======================================================================
// main

int main(int argc, char** argv)
{
  // psc_init(argc, argv);
  // FIXME restore whatever previous functionality there was with options
  int temp = 1;
  psc_init(temp, argv);

  run(argc, argv);

  psc_finalize();
  return 0;
}
