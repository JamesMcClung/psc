
#include "gtest/gtest.h"

#include <dim.hxx>
#include <psc/gtensor.h>
#include <psc/deposit.hxx>

#include <gtensor/reductions.h>
#include <mpi.h>

#include <limits>

using namespace gt::placeholders;

// ---------------------------------------------------------------------------

template <typename T>
struct DepositTest : ::testing::Test
{
  using real_t = typename T::real_t;
  using real3_t = Vec3<real_t>;
  using dim_t = typename T::dim_t;

  DepositTest() : ldims_{4, 4, 4}
  {
    if (dim_t::InvarX::value) {
      ldims_[0] = 1;
    }
  }

  template <typename F>
  void test_charge(real3_t x, const F& rho_ref)
  {
    const real_t eps = std::numeric_limits<real_t>::epsilon();

    auto flds = gt::zeros<real_t>(ldims_);
    psc::DepositNc<real_t, dim_t> deposit;
    deposit(flds, x);

    EXPECT_LT(gt::norm_linf(flds - rho_ref), eps) << flds;
  }

  gt::shape_type<3> ldims_;
};

template <typename R, typename D>
class DepositTestConfig
{
public:
  using real_t = R;
  using dim_t = D;
};

using DepositTestTypes = ::testing::Types<DepositTestConfig<double, dim_yz>,
                                          DepositTestConfig<double, dim_xyz>>;

TYPED_TEST_SUITE(DepositTest, DepositTestTypes);

TYPED_TEST(DepositTest, ChargeCenter)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;
  using dim_t = typename self_type::dim_t;

  real3_t x = {1.5, 1.5, 1.5};
  auto rho_ref = gt::zeros<real_t>(this->ldims_);
  if (std::is_same<dim_t, dim_xyz>::value) {
    // clang-format off
    rho_ref.view(1, _all, _all) = gt::gtensor<real_t, 2>({{0., 0.   , 0.   , 0.},
                                                          {0., 0.125, 0.125, 0.},
                                                          {0., 0.125, 0.125, 0.},
                                                          {0., 0.   , 0.   , 0.}});
    rho_ref.view(2, _all, _all) = gt::gtensor<real_t, 2>({{0., 0.   , 0.   , 0.},
                                                          {0., 0.125, 0.125, 0.},
                                                          {0., 0.125, 0.125, 0.},
                                                          {0., 0.   , 0.   , 0.}});
    // clang-format on
  } else if (std::is_same<dim_t, dim_yz>::value) {
    // clang-format off
    rho_ref.view(0, _all, _all) = gt::gtensor<real_t, 2>({{0., 0.  , 0.  , 0.},
                                                          {0., 0.25, 0.25, 0.},
                                                          {0., 0.25, 0.25, 0.},
                                                          {0., 0.  , 0.  , 0.}});
    // clang-format on
  }
  this->test_charge(x, rho_ref);
}

TYPED_TEST(DepositTest, ChargeLowerLeft)
{
  using self_type = DepositTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using real3_t = typename self_type::real3_t;
  using dim_t = typename self_type::dim_t;

  real3_t x = {1., 1., 1.};
  auto rho_ref = gt::zeros<real_t>(this->ldims_);
  if (std::is_same<dim_t, dim_xyz>::value) {
    // clang-format off
    rho_ref.view(1, _all, _all) = gt::gtensor<real_t, 2>({{0., 0., 0., 0.},
                                                          {0., 1., 0., 0.},
                                                          {0., 0., 0., 0.},
                                                          {0., 0., 0., 0.}});
    // clang-format on
  } else if (std::is_same<dim_t, dim_yz>::value) {
    // clang-format off
    rho_ref.view(0, _all, _all) = gt::gtensor<real_t, 2>({{0., 0., 0., 0.},
                                                          {0., 1., 0., 0.},
                                                          {0., 0., 0., 0.},
                                                          {0., 0., 0., 0.}});
    // clang-format on
  }
  this->test_charge(x, rho_ref);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
