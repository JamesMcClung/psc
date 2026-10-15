#include "neumann_test_util.hxx"

// v = 2 / sqrt(5) ~ .894 and dt ~ .53, so particles move ~.47 cells per step

TEST(NeumannBoundaryInjectorTest, InwardsLo)
{
  auto ys = run<LoHi::Lo>({{.75, 2.}}, 1);
  ASSERT_EQ(ys.size(), 3);
  EXPECT_NEAR(ys[0], ys[2] - 1., 1e-10); // copy, one cell behind the ion
  EXPECT_EQ(ys[1], .75);                 // electron
  EXPECT_GT(ys[2], 1.);                  // ion
}

TEST(NeumannBoundaryInjectorTest, InwardsHi)
{
  auto ys = run<LoHi::Hi>({{7.25, -2.}}, 1);
  ASSERT_EQ(ys.size(), 3);
  EXPECT_LT(ys[0], 7.);                  // ion
  EXPECT_EQ(ys[1], 7.25);                // electron
  EXPECT_NEAR(ys[2], ys[0] + 1., 1e-10); // copy, one cell behind the ion
}

// note: each run creates a psc integrator, and creating more than 4 in one
// process currently crashes at exit, so keep the number of tests small

TEST(NeumannBoundaryInjectorTest, NoCopies)
{
  auto ys = run<LoHi::Lo>({{.25, 1.},    // stays in edge cell
                           {.25, -2.},   // leaves domain
                           {7.25, -2.}}, // moves inwards at other boundary
                          1);
  ASSERT_EQ(ys.size(), 5); // outgoing ion was dropped
}

TEST(NeumannBoundaryInjectorTest, ManySteps)
{
  // each copy is itself copied when it leaves the edge cell; resampling with a
  // high temperature checks that the copies still conserve charge
  std::unique_ptr<Grid_t> grid{setupGrid(.75)};
  auto ys = run<LoHi::Lo>({{.75, 2.}}, 6,
                          NeumannResamplerMaxwellian{grid->kinds, {9., 9.}});
  ASSERT_GT(ys.size(), 3);
}

TEST(NeumannResamplerTest, Rotate)
{
  NeumannResamplerRotate resampler;
  psc::particle::Inject prt{{0, 0, 0}, {.3, .5, .4}, 1, KIND_ION};
  auto u = resampler.resample(prt, 1);
  EXPECT_EQ(u[1], .5);
  EXPECT_NEAR(sqr(u[0]) + sqr(u[2]), .25, 1e-12);
  EXPECT_NE(u[0], .3);
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
