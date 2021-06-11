#include "InsertionBuffer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace dyden {

TEST(InsertionBufferTest, ctor) {
  InsertionBuffer<int> a;
}

TEST(InsertionBufferTest, addValue) {
  InsertionBuffer<int> uut;

  for (int i = 0; i < 20; i++) {
    uut.addValue(i);
  }
}

TEST(InsertionBufferTest, flush) {
  InsertionBuffer<int> uut(25);

  int total_flushed = 0;
  int next_val = 1;

  for (int i = 0; i < 20; i++) {
    EXPECT_EQ(uut.addValue(next_val++), i + 1);
  }

  for (auto it = uut.lockedIterator(); it; ++it) {
    total_flushed += *it;
  }

  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(uut.addValue(next_val++), i + 1);
  }

  for (auto it = uut.lockedIterator(); it; ++it) {
    total_flushed += *it;
  }

  EXPECT_EQ(total_flushed, next_val * (next_val - 1) / 2);
}

TEST(InsertionBufferTest, overflow) {
  InsertionBuffer<int> uut(/*buffer_size=*/10);

  for (int i = 0; i < 20; i++) {
    uut.addValue(i);
  }

  int num_flushed = 0;
  int total_flushed = 0;
  for (auto it = uut.lockedIterator(); it; ++it) {
    num_flushed++;
    total_flushed += *it;
  }
  EXPECT_EQ(num_flushed, 9);
  // Only the last 9 values should be in the total, so subtract the first 11
  // from the sum.
  EXPECT_EQ(total_flushed, 20 * 19 / 2 - 11 * 10 / 2);
}

}  // namespace dhist
