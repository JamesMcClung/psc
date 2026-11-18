
#include <gtest/gtest.h>

#include "psc_fields_c.h"
#include "fields.hxx"
#include "interpolate.hxx"

#include <memory>

using real_t = MfieldsC::real_t;

namespace
{

// distinct linear profiles per component, f_m(y) = a_m + b_m y
double a(int m) { return 1. + m; }
double b(int m) { return .25 * (m + 1) * (m % 2 ? -1. : 1.); }

// offset of component m's y location from its index, in cells, on the Yee
// grid: nodes (0) for EX, EZ, HY; cell centers (.5) for EY, HX, HZ
double y_offset(int m) { return (m == EY || m == HX || m == HZ) ? .5 : 0.; }

double profile(int m, double y) { return a(m) + b(m) * y; }

struct InterpolateDimYTest : ::testing::Test
{
  InterpolateDimYTest()
  {
    Int3 gdims{1, 8, 1};
    auto domain = Grid_t::Domain{gdims, Double3(gdims)}; // dx = 1
    auto norm = Grid_t::Normalization{};
    grid_.reset(new Grid_t{domain, {}, {}, norm, 1., -1, {0, 2, 0}});
    mflds_.reset(new MfieldsC{*grid_, NR_FIELDS, grid_->ibn});

    for (int m = EX; m <= HZ; m++) {
      for (int jy = -2; jy < 10; jy++) {
        (*mflds_)(m, 0, jy, 0, 0) = profile(m, jy + y_offset(m));
      }
    }
  }

  template <typename OPT_IP>
  void check(double y, bool ec)
  {
    auto EM = make_Fields3d<dim_y>((*mflds_)[0]);
    using IP = InterpolateEM<decltype(EM), OPT_IP, dim_y>;
    IP ip;
    real_t xm[3] = {.5, real_t(y), .5};
    ip.set_coeffs(xm);

    // with EC, EY, HX, HZ aren't interpolated along y but taken from the cell
    int l = ip.cy.g.l;
    auto expected = [&](int m) {
      bool from_cell = ec && (m == EY || m == HX || m == HZ);
      return from_cell ? profile(m, l + y_offset(m)) : profile(m, y);
    };

    const double eps = 1e-12;
    EXPECT_NEAR(ip.ex(EM), expected(EX), eps) << "y " << y;
    EXPECT_NEAR(ip.ey(EM), expected(EY), eps) << "y " << y;
    EXPECT_NEAR(ip.ez(EM), expected(EZ), eps) << "y " << y;
    EXPECT_NEAR(ip.hx(EM), expected(HX), eps) << "y " << y;
    EXPECT_NEAR(ip.hy(EM), expected(HY), eps) << "y " << y;
    EXPECT_NEAR(ip.hz(EM), expected(HZ), eps) << "y " << y;
  }

  std::unique_ptr<Grid_t> grid_;
  std::unique_ptr<MfieldsC> mflds_;
};

const double test_ys[] = {2.0, 2.3, 2.5, 2.9};

} // namespace

TEST_F(InterpolateDimYTest, Std1stDimY)
{
  for (double y : test_ys) {
    check<opt_ip_1st>(y, false);
  }
}

TEST_F(InterpolateDimYTest, Ec1stDimY)
{
  for (double y : test_ys) {
    check<opt_ip_1st_ec>(y, true);
  }
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
