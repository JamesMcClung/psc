#include "neumann_test_util.hxx"

// Copies made deeper than the edge cell (copy_depth > 1). A separate binary
// from test_neumann_boundary_injector because one binary may create at most 4
// integrators.

// v = 2 / sqrt(5) ~ .894 and dt ~ .53, so particles move ~.47 cells per step

TEST(NeumannBoundaryInjectorDepthTest, Lo)
{
  // staggered so that a copy of the first ion would sit elsewhere
  auto ys = run<LoHi::Lo>({{.6, 2.},    // crosses y = 1: not copied
                           {2.75, 2.}}, // crosses y = 3: copied
                          1, {}, 3);
  ASSERT_EQ(ys.size(), 5);
  EXPECT_LT(ys[0], 1.);                  // copy, three cells behind the ion
  EXPECT_EQ(ys[1], .6);                  // electron
  EXPECT_GT(ys[2], 1.);                  // first ion
  EXPECT_EQ(ys[3], 2.75);                // electron
  EXPECT_NEAR(ys[0], ys[4] - 3., 1e-10); // second ion
}

TEST(NeumannBoundaryInjectorDepthTest, Hi)
{
  auto ys = run<LoHi::Hi>({{5.25, -2.}}, 1, {}, 3); // crosses y = 5: copied
  ASSERT_EQ(ys.size(), 3);
  EXPECT_LT(ys[0], 5.);                  // ion
  EXPECT_EQ(ys[1], 5.25);                // electron
  EXPECT_NEAR(ys[2], ys[0] + 3., 1e-10); // copy, three cells behind the ion
}

TEST(NeumannBoundaryInjectorDepthTest, ManySteps)
{
  // copies of copies, resampled hot, still conserve charge
  std::unique_ptr<Grid_t> grid{setupGrid(.75)};
  auto ys = run<LoHi::Lo>({{2.75, 2.}}, 10,
                          NeumannResamplerMaxwellian{grid->kinds, {9., 9.}}, 3);
  ASSERT_GT(ys.size(), 3);
}

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int rc = RUN_ALL_TESTS();
  MPI_Finalize();
  return rc;
}
