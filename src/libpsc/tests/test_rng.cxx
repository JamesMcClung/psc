#include <gtest/gtest.h>

#include <set>

#include "rng.hxx"

TEST(Rng, ProcessSeedWithoutBaseIsHistorical)
{
  EXPECT_EQ(rng::detail::process_seed(0, std::nullopt), 10000);
  EXPECT_EQ(rng::detail::process_seed(3, std::nullopt), 40000);
}

TEST(Rng, ProcessSeedWithBaseIsDeterministic)
{
  EXPECT_EQ(rng::detail::process_seed(2, 7), rng::detail::process_seed(2, 7));
}

TEST(Rng, ProcessSeedWithBaseVariesWithBaseAndRank)
{
  std::set<int> seeds;
  for (int base = 1; base <= 10; base++) {
    for (int rank = 0; rank < 8; rank++) {
      int seed = rng::detail::process_seed(rank, base);
      EXPECT_GE(seed, 0);
      seeds.insert(seed);
    }
  }
  EXPECT_EQ(seeds.size(), 80);
}

// creates the process generator, so it must stay the last test
TEST(Rng, SetBaseSeedAfterFirstDrawThrows)
{
  rng::Uniform<double> uniform;
  uniform.get();
  EXPECT_THROW(rng::set_base_seed(1), std::logic_error);
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
