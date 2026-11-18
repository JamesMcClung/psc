#include "gtest/gtest.h"

#include <psc.hxx>

#include "../psc_config.hxx"
#include "../psc_bnd_fields/conducting_wall.hxx"

using Dim = dim_yz;
using PscConfig = PscConfig1vbecDouble<Dim>;
using MfieldsState = PscConfig::MfieldsState;

using psc::Axis;
using psc::bnd::LoHi;
using psc::bnd::field::ConductingWall;

namespace
{

const int n_ghosts = 2;

Grid_t makeGrid()
{
  auto domain =
    Grid_t::Domain{{1, 8, 6}, {1.0, 8.0, 6.0}, {0, 0, 0}, {1, 1, 1}};
  auto bc = psc::grid::BC{
    {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_CONDUCTING_WALL},
    {BND_FLD_PERIODIC, BND_FLD_CONDUCTING_WALL, BND_FLD_CONDUCTING_WALL},
    {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_REFLECTING},
    {BND_PRT_PERIODIC, BND_PRT_REFLECTING, BND_PRT_REFLECTING}};
  auto kinds = Grid_t::Kinds{};
  auto norm = Grid_t::Normalization{};
  return Grid_t{
    domain, bc, kinds, norm, 1e-2, -1, n_ghosts * Dim::get_noninvariant_mask()};
}

// fills every E and H value, ghosts included, with a distinct nonzero number
void fill(MfieldsState& mflds)
{
  Int3 ib = mflds.ib(), ie = mflds.ib() + mflds.im();
  for (int m = EX; m <= HZ; m++) {
    for (Int3 i3 : VecRange(ib, ie)) {
      mflds(m, i3[0], i3[1], i3[2], 0) = 1 + 1000 * m + 100 * (i3[0] - ib[0]) +
                                         10 * (i3[1] - ib[1]) + (i3[2] - ib[2]);
    }
  }
}

// Checks the wall's ghosts along axis `d`. The wall is at node `w` (0 or
// ldims[d]), and `out` = -1 or +1 points out of the domain. Node-centered
// values at w + out*j mirror w - out*j. Cell-centered values (index i is at
// i + 1/2) at the j-th cell outside the wall mirror the j-th cell inside it.
void check_wall(MfieldsState& mflds, int d, LoHi lohi)
{
  const Grid_t& grid = mflds.grid();
  int w = lohi == LoHi::Lo ? 0 : grid.ldims[d];
  int out = lohi == LoHi::Lo ? -1 : 1;
  // on the hi side, node-centered values have one less ghost
  int n_nc = lohi == LoHi::Lo ? n_ghosts : n_ghosts - 1;

  Axis ax = d;
  int e_n = EX + d, e_t1 = EX + ax.next(), e_t2 = EX + ax.prev();
  int h_n = HX + d, h_t1 = HX + ax.next(), h_t2 = HX + ax.prev();

  // index of the j-th node or cell outside (j > 0) or inside (j <= 0)
  auto nc = [&](int j) { return w + out * j; };
  auto cc = [&](int j) { return lohi == LoHi::Lo ? w - j : w + j - 1; };

  Int3 ib = mflds.ib(), ie = mflds.ib() + mflds.im();
  ib[d] = 0;
  ie[d] = 1;
  for (Int3 i3 : VecRange(ib, ie)) {
    auto F = [&](int m, int i) -> double {
      Int3 idx = i3;
      idx[d] = i;
      return mflds(m, idx[0], idx[1], idx[2], 0);
    };

    for (int m : {e_t1, e_t2}) {
      EXPECT_EQ(F(m, nc(0)), 0.0) << "m " << m << " i3 " << i3;
      for (int j = 1; j <= n_nc; j++) {
        EXPECT_EQ(F(m, nc(j)), F(m, nc(-j))) << "m " << m << " j " << j;
      }
    }
    for (int j = 1; j <= n_nc; j++) {
      EXPECT_EQ(F(h_n, nc(j)), F(h_n, nc(-j))) << "j " << j;
    }
    for (int j = 1; j <= n_ghosts; j++) {
      EXPECT_EQ(F(e_n, cc(j)), -F(e_n, cc(1 - j))) << "j " << j;
      for (int m : {h_t1, h_t2}) {
        EXPECT_EQ(F(m, cc(j)), -F(m, cc(1 - j))) << "m " << m << " j " << j;
      }
    }
  }
}

void test_wall(int d, LoHi lohi)
{
  auto grid = makeGrid();
  MfieldsState mflds{grid};
  fill(mflds);

  ConductingWall<Dim, MfieldsState> wall{d, lohi};
  wall.apply_e_bcs(mflds);
  wall.apply_h_bcs(mflds);

  check_wall(mflds, d, lohi);
}

} // namespace

TEST(ConductingWallTest, LoY) { test_wall(1, LoHi::Lo); }
TEST(ConductingWallTest, HiY) { test_wall(1, LoHi::Hi); }
TEST(ConductingWallTest, LoZ) { test_wall(2, LoHi::Lo); }
TEST(ConductingWallTest, HiZ) { test_wall(2, LoHi::Hi); }

int main(int argc, char** argv)
{
  psc_init(argc, argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  psc_finalize();
  return rc;
}
