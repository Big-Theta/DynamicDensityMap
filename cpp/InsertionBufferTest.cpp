#include "cpp/InsertionBuffer.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace dhist {

TEST(InsertionBufferTest, ctor) {
  InsertionBuffer</*T=*/int, /*kThreadsafe=*/false> a;
  InsertionBuffer</*T=*/int, /*kThreadsafe=*/true> b;
  EXPECT_TRUE(true);
}

TEST(InsertionBufferTest, addValue) {
  InsertionBuffer</*T=*/int, /*kThreadsafe=*/true> uut;

  for (int i = 0; i < 20; i++) {
    uut.addValue(i);
  }
}

TEST(InsertionBufferTest, flush) {
  InsertionBuffer</*T=*/int, /*kThreadsafe=*/true> uut(25);

  int total_flushed = 0;
  int next_val = 1;

  for (int i = 0; i < 20; i++) {
    uut.addValue(next_val++);
  }

  for (auto it = uut.begin(); it != uut.end(); it++) {
    total_flushed += *it;
  }

  for (int i = 0; i < 10; i++) {
    uut.addValue(next_val++);
  }

  for (auto it = uut.begin(); it != uut.end(); it++) {
    total_flushed += *it;
  }

  EXPECT_EQ(total_flushed, next_val * (next_val - 1) / 2);
}

}  // namespace dhist
