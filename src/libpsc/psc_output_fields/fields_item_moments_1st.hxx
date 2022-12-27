
#pragma once

#include <psc/moment.hxx>
#include "fields_item.hxx"

#include <cmath>

template <typename S, typename D>
using Moment_n_1st =
  ItemMoment<psc::moment::moment_n<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moment_v_1st =
  ItemMoment<psc::moment::moment_v<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moment_p_1st =
  ItemMoment<psc::moment::moment_p<psc::deposit::code::Deposit1stCc, D>, S>;

template <typename S, typename D>
using Moment_T_1st =
  ItemMoment<psc::moment::moment_T<psc::deposit::code::Deposit1stCc, D>, S>;

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename MP, typename S, typename D>
class Moments_1st : public ItemMomentCRTP<Moments_1st<MP, S, D>, S>
{
public:
  using Base = ItemMomentCRTP<Moments_1st<MP, S, D>, S>;
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, D>;

  using Base::Base;
};

#ifdef USE_CUDA

#include "../libpsc/cuda/mparticles_cuda.hxx"
#include "psc_particles_single.h"

// ======================================================================
// Moments_1st
//
// all moments calculated at once
// FIXME: add KE

template <typename BS, typename D>
class Moments_1st<MparticlesCuda<BS>, MfieldsSingle::Storage, D>
  : public ItemMomentCRTP<
      Moments_1st<MparticlesCuda<BS>, MfieldsSingle::Storage, D>,
      MfieldsSingle::Storage>
{
public:
  using Base =
    ItemMomentCRTP<Moments_1st<MparticlesCuda<BS>, MfieldsSingle::Storage, D>,
                   MfieldsSingle::Storage>;
  using dim_t = D;
  using Mparticles = MparticlesCuda<BS>;
  using Mfields = MfieldsSingle;
  using real_t = Mfields::real_t;
  using moment_type =
    psc::moment::moment_all<psc::deposit::code::Deposit1stCc, dim_t>;

  explicit Moments_1st(const Grid_t& grid) : Base{grid} {}

  auto operator()(const Mparticles& _mprts)
  {
    static int pr, pr_A, pr_B, pr_C;
    if (!pr) {
      pr = prof_register("Moments_1st cuda", 1., 0, 0);
      pr_A = prof_register("Moments_1st get", 1., 0, 0);
      pr_B = prof_register("Moments_1st process", 1., 0, 0);
      pr_C = prof_register("Moments_1st addg", 1., 0, 0);
    }

    prof_start(pr);
    prof_start(pr_A);
    auto& mprts = const_cast<Mparticles&>(_mprts);
    auto&& h_mprts = mprts.template get_as<MparticlesSingle>();
    prof_stop(pr_A);

    prof_start(pr_B);
    moment_type()(Base::mres_gt_, Base::mres_ib_, h_mprts);
    prof_stop(pr_B);

    prof_start(pr_C);
    Base::bnd_.add_ghosts(mprts.grid(), Base::mres_gt_, Base::mres_ib_);
    prof_stop(pr_C);

    mprts.put_as(h_mprts, MP_DONT_COPY);
    prof_stop(pr);
    return Base::mres_gt_;
  }
};

#endif
