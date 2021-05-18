// MIT License
//
// Copyright (c) 2020 Logan Evans
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include <google/protobuf/util/time_util.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

#include "DensityMap.pb.h"

namespace dhist {

class Kernel {
  
};

template <bool kUseDecay = true, bool kThreadsafe = true>
class DynamicKDE {
 public:
  DynamicKDE(size_t num_kernels, double decay_rate = 0.0,
             int refresh_interval = 512) {}

 private:
  const double decay_rate_;
  const int refresh_interval_;
  uint64_t generation_;
  uint64_t refresh_generation_;
  double total_count_;

  int insert_queue_begin_;
  int insert_queue_end_;
  int insert_queue_to_flush_;

  std::vector<Kernel> kernels_;
};

}  // namespace dhist
