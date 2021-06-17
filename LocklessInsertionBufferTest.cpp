#include "LocklessInsertionBuffer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace dyden {

TEST(LocklessInsertionBufferTest, ctor) { LocklessInsertionBuffer<int> a; }

TEST(LocklessInsertionBufferTest, addValue) {
  LocklessInsertionBuffer<int> uut;

  for (int i = 0; i < 20; i++) {
    uut.addValue(i);
  }
}

TEST(LocklessInsertionBufferTest, flush) {
  LocklessInsertionBuffer<int> uut(32);

  int total_flushed = 0;
  int next_val = 1;

  for (int i = 0; i < 30; i++) {
    uut.addValue(next_val++);
  }

  for (auto it = uut.lockedIterator(); it; ++it) {
    total_flushed += *it;
  }

  EXPECT_EQ(total_flushed, next_val * (next_val - 1) / 2);

  for (int i = 0; i < 10; i++) {
    uut.addValue(next_val++);
  }

  for (auto it = uut.lockedIterator(); it; ++it) {
    total_flushed += *it;
  }

  EXPECT_EQ(total_flushed, next_val * (next_val - 1) / 2);
}

TEST(LocklessInsertionBufferTest, overflow) {
  LocklessInsertionBuffer<int> uut(/*buffer_size=*/16);

  for (int i = 0; i < 20; i++) {
    uut.addValue(i);
  }

  int num_flushed = 0;
  int total_flushed = 0;
  for (auto it = uut.lockedIterator(); it; ++it) {
    num_flushed++;
    total_flushed += *it;
  }
  EXPECT_EQ(num_flushed, 16);
  EXPECT_EQ(total_flushed, 20 * 19 / 2 - 4 * 3 / 2);
}

}  // namespace dyden
