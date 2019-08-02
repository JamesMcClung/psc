
#pragma once

#include "fields_item.hxx"
#include "bnd_cuda_3_impl.hxx"
#include "psc_fields_cuda.h"
#include "cuda_moments.cuh"

template<typename BS>
struct cuda_mparticles;

// ======================================================================
// Moment_rho_1st_nc_cuda

template<typename _Mparticles, typename dim>
struct Moment_rho_1st_nc_cuda
{
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda3<Mfields>;
  
  constexpr static const char* name = "rho_1st_nc";
  static int n_comps(const Grid_t&) { return 1; }
  static std::vector<std::string> fld_names() { return { "rho_nc_cuda" }; } // FIXME
  constexpr static int flags = 0;

  Moment_rho_1st_nc_cuda(const Grid_t& grid)
    : mres_{grid, n_comps(grid), grid.ibn},
      bnd_{grid, grid.ibn}
  {}

  void operator()(Mparticles& mprts)
  {
    auto& cmprts = *mprts.cmprts();
    cuda_mfields *cmres = mres_.cmflds();
    
    mres_.zero();
    CudaMoments1stNcRho<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, cmres);
    bnd_.add_ghosts(mres_, 0, mres_.n_comps());
  }

  Mfields& result() { return mres_; }

private:
  Mfields mres_;
  Bnd bnd_;
};

// ======================================================================
// Moment_n_1st_cuda

template<typename _Mparticles, typename dim>
struct Moment_n_1st_cuda
{
  using Mparticles = _Mparticles;
  using Mfields = MfieldsCuda;
  using Bnd = BndCuda3<Mfields>;
  
  constexpr static const char* name = "n_1st";
  static int n_comps(const Grid_t&) { return 1; }
  static std::vector<std::string> fld_names() { return { "n_1st_cuda" }; }
  constexpr static int flags = POFI_BY_KIND;

  Moment_n_1st_cuda(const Grid_t& grid)
    : mres_{grid, n_comps(grid), grid.ibn},
      bnd_{grid, grid.ibn}
  {}

  void operator()(Mparticles& mprts)
  {
    auto& cmprts = *mprts.cmprts();
    cuda_mfields *cmres = mres_.cmflds();
    
    mres_.zero();
    CudaMoments1stNcN<cuda_mparticles<typename Mparticles::BS>, dim> cmoments;
    cmoments(cmprts, cmres);
    bnd_.add_ghosts(mres_, 0, mres_.n_comps());
  }

  Mfields& result() { return mres_; }

private:
  Mfields mres_;
  Bnd bnd_;
};

