

#include "gtest/gtest.h"

#include "testing.hxx"

#include <psc/deposit.hxx>
#include "grid.hxx"
#include "fields.hxx"

#include <limits>

// ---------------------------------------------------------------------------

template <typename D>
class CalcDivNc
{
public:
  using dim_t = D;

  template <typename F>
  auto operator()(const F& flds, dim_y tag)
  {
    return gt::eval(flds.view(_all, _s(1, _), _all, 1) -
                    flds.view(_all, _s(0, -1), _all, 1));
  }

  template <typename F>
  auto operator()(const F& flds, dim_yz tag)
  {
    return gt::eval((flds.view(_all, _s(1, _), _s(1, _), 1) -
                     flds.view(_all, _s(0, -1), _s(1, _), 1)) +
                    (flds.view(_all, _s(1, _), _s(1, _), 2) -
                     flds.view(_all, _s(1, _), _s(0, -1), 2)));
  }

  template <typename F>
  auto operator()(const F& flds, dim_xyz tag)
  {
    return gt::eval((flds.view(_s(1, _), _s(1, _), _s(1, _), 0) -
                     flds.view(_s(0, -1), _s(1, _), _s(1, _), 0)) +
                    (flds.view(_s(1, _), _s(1, _), _s(1, _), 1) -
                     flds.view(_s(1, _), _s(0, -1), _s(1, _), 1)) +
                    (flds.view(_s(1, _), _s(1, _), _s(1, _), 2) -
                     flds.view(_s(1, _), _s(1, _), _s(0, -1), 2)));
  }

  template <typename F>
  auto operator()(const F& flds)
  {
    return (*this)(flds, dim_t{});
  }
};

// ---------------------------------------------------------------------------

template <typename T>
struct CurrentDepositionTest : ::testing::Test
{
  using real_t = typename T::real_t;
  using Real3 = Vec3<real_t>;
  using dim_t = typename T::dim_t;
  using Current = typename T::Current;

  CurrentDepositionTest()
  {
    double dt = 1.;
    Int3 gdims = {4, 4, 4};
    if (dim_t::InvarX::value) {
      gdims[0] = 1;
    }
    if (dim_t::InvarZ::value) {
      gdims[2] = 1;
    }
    Double3 length;
    for (int d = 0; d < 3; d++) {
      length[d] = double(gdims[d]);
    } // make length == gdims, so that dx = 1
    psc::grid::Domain<double> grid_domain{gdims, length};
    psc::grid::BC grid_bc;
    Grid_t::Normalization grid_norm;
    std::vector<Grid_t::Kind> kinds;

    grid_.reset(
      new Grid_t{grid_domain, grid_bc, kinds, grid_norm, dt, 1, this->ibn});

    ldims_ = Int3(grid().ldims);
  }

  const Grid_t& grid() const { return *grid_; }

  // lg is the cell the particle starts in, as the pusher passes it; by default
  // the one containing xm
  auto calc_current(Real3 xm, Real3 xp, Real3 vxi, Int3 lg = Int3{-1, -1, -1})
  {
    if (lg == Int3{-1, -1, -1}) {
      lg = xm.fint();
    }
    using fields_view_t = typename T::fields_view_t;
    using fields_t = typename T::fields_t;

    Current curr{grid()};

    auto flds =
      gt::zeros<real_t>(gt::shape(ldims_[0], ldims_[1], ldims_[2], NR_FIELDS));

    // FIXME, to_kernel() is a hack to get a gtensor_view
    auto flds_view = fields_view_t({0, 0, 0}, flds.to_kernel());
    fields_t J{flds_view};

    real_t qni_wni = 1.;
    Int3 lf = xp.fint();
    curr.calc_j(J, xm, xp, lf, lg, qni_wni, vxi);
    return flds;
  }

  template <typename F>
  void check_continuity(Real3 xm, Real3 xp, const F& flds)
  {
    const real_t eps = 2. * std::numeric_limits<real_t>::epsilon();

    Int3 ib{0, 0, 0};
    auto rho_m = gt::zeros<real_t>(ldims_);
    auto rho_p = gt::zeros<real_t>(ldims_);
    psc::deposit::norm::nc<dim_t>(rho_m, ib, xm, 1.);
    psc::deposit::norm::nc<dim_t>(rho_p, ib, xp, 1.);
    gt::gtensor<real_t, 3> d_rho;
    if (std::is_same<dim_t, dim_y>::value) {
      d_rho = (rho_p - rho_m).view(_all, _s(1, _), _all);
    } else if (std::is_same<dim_t, dim_yz>::value) {
      d_rho = (rho_p - rho_m).view(_all, _s(1, _), _s(1, _));
    } else {
      d_rho = (rho_p - rho_m).view(_s(1, _), _s(1, _), _s(1, _));
    };

    CalcDivNc<dim_t> div;
    auto div_j = div(flds.view(_all, _all, _all, _s(JXI, JXI + 3)));
    // std::cout << "d_rho\n" << d_rho.view(0, _all, _all) << "\n";
    // std::cout << "div_j\n" << div_j.view(0, _all, _all) << "\n";
    EXPECT_LT(gt::norm_linf(d_rho + div_j), eps);
  }

  template <typename F>
  void test_current(Real3 xm, Real3 xp, Real3 vxi, const F& jxi_ref,
                    const F& jyi_ref, const F& jzi_ref)
  {
    if (!std::is_same<dim_t, dim_yz>::value) {
      return;
    }
    auto flds = calc_current(xm, xp, vxi);
    check_continuity(xm, xp, flds);

    const real_t eps = 2. * std::numeric_limits<real_t>::epsilon();
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all, JXI) - jxi_ref), eps)
      << flds.view(0, _all, _all, JXI);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all, JYI) - jyi_ref), eps)
      << flds.view(0, _all, _all, JYI);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, _all, JZI) - jzi_ref), eps)
      << flds.view(0, _all, _all, JZI);
  }

  // 1d (dim_y) refs, indexed by y
  template <typename F>
  void test_current_1d(Real3 xm, Real3 xp, Real3 vxi, const F& jxi_ref,
                       const F& jyi_ref, const F& jzi_ref)
  {
    if (!std::is_same<dim_t, dim_y>::value) {
      return;
    }
    auto flds = calc_current(xm, xp, vxi);
    check_continuity(xm, xp, flds);

    const real_t eps = 2. * std::numeric_limits<real_t>::epsilon();
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, 0, JXI) - jxi_ref), eps)
      << flds.view(0, _all, 0, JXI);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, 0, JYI) - jyi_ref), eps)
      << flds.view(0, _all, 0, JYI);
    EXPECT_LT(gt::norm_linf(flds.view(0, _all, 0, JZI) - jzi_ref), eps)
      << flds.view(0, _all, 0, JZI);
  }

  void test_current(Real3 xm, Real3 xp, Real3 vxi, Int3 lg = Int3{-1, -1, -1})
  {
    auto flds = calc_current(xm, xp, vxi, lg);
    check_continuity(xm, xp, flds);
  }

  Int3 ibn = {0, 0, 0};
  std::unique_ptr<Grid_t> grid_;
  Int3 ldims_;
};

template <typename R, typename D,
          template <typename, typename, typename> class C>
class CurrentDepositionTestConfig
{
public:
  using real_t = R;
  using dim_t = D;
  using fields_view_t = SArrayView<real_t, gt::space::host>;
  using fields_t = curr_cache_t<fields_view_t>;
  using Current = C<opt_order_1st, dim_t, fields_t>;
};

using CurrentDepositionTestTypes =
  ::testing::Types< // CurrentDepositionTestConfig<float, dim_yz,
                    // Current1vbSplit>,
    CurrentDepositionTestConfig<double, dim_y, Current1vbSplit>,
    CurrentDepositionTestConfig<double, dim_yz, Current1vbSplit>,
    CurrentDepositionTestConfig<double, dim_yz, CurrentZigzag>,
    CurrentDepositionTestConfig<double, dim_xyz, CurrentZigzag>,
    CurrentDepositionTestConfig<double, dim_xyz, Current1vbSplit>>;

TYPED_TEST_SUITE(CurrentDepositionTest, CurrentDepositionTestTypes);

TYPED_TEST(CurrentDepositionTest, CurrentNotMoving)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1., 1.}, xp = {.5, 1., 1.};
  Real3 vxi = {0., 0., 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentY)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1., 1.}, xp = {.5, 1.2, 1.};
  Real3 vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0. , 0., 0.},
                                    {0., 0.2, 0., 0.},
                                    {0., 0. , 0., 0.},
                                    {0., 0. , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentYShift)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1., 1.7}, xp = {.5, 1.2, 1.7};
  Real3 vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0., 0.},
                                    {0., 0.06, 0., 0.},
                                    {0., 0.14, 0., 0.},
                                    {0., 0.  , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentZ)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1., 1.3}, xp = {.5, 1., 1.6};
  Real3 vxi = {0., 0., 3.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0.3, 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentYCross)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1.9, 1.}, xp = {.5, 2.1, 1.};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  Real3 vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0. , 0., 0.},
                                    {0., 0.1, 0.1, 0.},
                                    {0., 0. , 0., 0.},
                                    {0., 0. , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentYCrossShift)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1.9, 1.3}, xp = {.5, 2.1, 1.3};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  Real3 vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.07, 0.07, 0.},
                                    {0., 0.03, 0.03, 0.},
                                    {0., 0.  , 0.  , 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentYZ)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1.2, 1.1}, xp = {.5, 1.4, 1.4};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  Real3 vxi = {0., .2, .3};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0., 0.},
                                    {0., 0.15, 0., 0.},
                                    {0., 0.05, 0., 0.},
                                    {0., 0.  , 0., 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.21, 0.09, 0.},
                                    {0., 0.  , 0.  , 0.},
                                    {0., 0.  , 0.  , 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentYZCrossShift)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1.9, 1.3}, xp = {.5, 2.1, 1.4};
  Int3 lf = {0, 1, 1}, lg = {0, 1, 1};
  Real3 vxi = {0., .2, 0.};
  // clang-format off
  gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0.  , 0.},
                                    {0., 0.0675, 0.0625, 0.},
                                    {0., 0.0325, 0.0375, 0.},
                                    {0., 0.  , 0.  , 0.}};
  gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                    {0., 0.0025, 0.095, 0.0025},
                                    {0., 0., 0., 0.},
                                    {0., 0., 0., 0.}};
  // clang-format on
  this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentYZCrossYZ)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1.9, 1.6}, xp = {.5, 2.1, 2.2};
  Real3 vxi = {0., .2, 0.};
  if (std::is_same<typename TypeParam::Current,
                   Current1vbSplit<opt_order_1st, dim_yz,
                                   typename TypeParam::fields_t>>::value) {
    // clang-format off
    gt::gtensor<real_t, 2> jxi_ref = {{0., 0., 0., 0.},
                                      {0., 0., 0., 0.},
                                      {0., 0., 0., 0.},
                                      {0., 0., 0., 0.}};
    gt::gtensor<real_t, 2> jyi_ref = {{0., 0.  , 0.  , 0.},
                                      {0., 0.025, 1./600., 0.},
                                      {0., 0.075, .09+1./600, 0.},
                                      {0., 0., 2./300., 0.}};
    gt::gtensor<real_t, 2> jzi_ref = {{0., 0., 0., 0.},
                                      {0., 0.015, 0.38+1/300., 1./600.},
                                      {0., 0., 0.18 + 2./300., 4./300.},
                                      {0., 0., 0., 0.}};
    // clang-format on
    this->test_current(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
  } else { // CurrentZigzag
    this->test_current(xm, xp, vxi);
  }
}

TYPED_TEST(CurrentDepositionTest, CurrentX)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {1.3, 1., 1.}, xp = {1.6, 1., 1.};
  Real3 vxi = {.3, 0., 0.};
  this->test_current(xm, xp, vxi);
}

TYPED_TEST(CurrentDepositionTest, CurrentXY)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {1.1, 1.2, 1.3}, xp = {1.4, 1.8, 1.3};
  Real3 vxi = {.3, .6, 0.};
  this->test_current(xm, xp, vxi);
}

TYPED_TEST(CurrentDepositionTest, CurrentXYZ)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {1.1, 1.2, 1.3}, xp = {1.4, 1.6, 1.8};
  Real3 vxi = {.3, .4, .5};
  this->test_current(xm, xp, vxi);
}

TYPED_TEST(CurrentDepositionTest, CurrentXYZCrossXYZ)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {1.9, 1.8, 1.7}, xp = {2.2, 2.4, 2.6};
  Real3 vxi = {.3, .6, .9};
  this->test_current(xm, xp, vxi);
}

TYPED_TEST(CurrentDepositionTest, CurrentY1d)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1., .5}, xp = {.5, 1.2, .5};
  Real3 vxi = {.1, .2, .3};
  gt::gtensor<real_t, 1> jxi_ref = {0., .09, .01, 0.};
  gt::gtensor<real_t, 1> jyi_ref = {0., .2, 0., 0.};
  gt::gtensor<real_t, 1> jzi_ref = {0., .27, .03, 0.};
  this->test_current_1d(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentY1dCross)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 1.8, .5}, xp = {.5, 2.3, .5};
  Real3 vxi = {.1, .5, .3};
  gt::gtensor<real_t, 1> jxi_ref = {0., .004, .087, .009};
  gt::gtensor<real_t, 1> jyi_ref = {0., .2, .3, 0.};
  gt::gtensor<real_t, 1> jzi_ref = {0., .012, .261, .027};
  this->test_current_1d(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

TYPED_TEST(CurrentDepositionTest, CurrentY1dCrossBack)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using real_t = typename self_type::real_t;
  using Real3 = typename self_type::Real3;

  Real3 xm = {.5, 2.3, .5}, xp = {.5, 1.8, .5};
  Real3 vxi = {.1, -.5, .3};
  gt::gtensor<real_t, 1> jxi_ref = {0., .004, .087, .009};
  gt::gtensor<real_t, 1> jyi_ref = {0., -.2, -.3, 0.};
  gt::gtensor<real_t, 1> jzi_ref = {0., .012, .261, .027};
  this->test_current_1d(xm, xp, vxi, jxi_ref, jyi_ref, jzi_ref);
}

// the upper boundary injector starts particles on the outer face of the upper
// ghost cell, i.e., on the upper face of the cell they start in
template <typename C>
struct is_current_1vb_split : std::false_type
{};

template <typename O, typename D, typename F>
struct is_current_1vb_split<Current1vbSplit<O, D, F>> : std::true_type
{};

TYPED_TEST(CurrentDepositionTest, CurrentFromUpperFaceY)
{
  using self_type = CurrentDepositionTest<TypeParam>;
  using Real3 = typename self_type::Real3;
  if (!is_current_1vb_split<typename self_type::Current>::value) {
    return;
  }

  Real3 xm = {1.5, 3., 1.5}, xp = {1.5, 1.6, 1.5};
  Real3 vxi = {0., -1.4, 0.};
  this->test_current(xm, xp, vxi, {1, 2, 1});
}

// records the highest y index a deposit adds to, even if it adds 0
struct MaxYRecordingCurr
{
  using real_t = double;
  using value_type = double;
  void add(int m, int i, int j, int k, double val) const
  {
    *max_j = std::max(*max_j, j);
  }
  int* max_j;
};

TEST(Current1vbSplitTest, EndOnUpperFaceY)
{
  // OpenBoundaryOutflow moves a particle in the upper ghost cell (here [2, 3))
  // to its outer face, y = 3. The deposit is in cell 2, so it reaches node 3,
  // but no further.
  auto domain = Grid_t::Domain{{1, 2, 1}, {1., 2., 1.}};
  auto norm =
    Grid_t::Normalization{Grid_t::NormalizationParams::dimensionless()};
  Grid_t grid{domain, {}, {}, norm, .1, -1, {0, 2, 0}};
  Current1vbSplit<opt_order_1st, dim_y, MaxYRecordingCurr> current{grid};

  int max_j = std::numeric_limits<int>::min();
  Vec3<double> xm = {.5, 2.4, .5}, xp = {.5, 3., .5}, vxi = {0., 1., 0.};
  current.calc_j(MaxYRecordingCurr{&max_j}, xm, xp, {}, {}, 1., vxi);
  EXPECT_EQ(max_j, 3);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
