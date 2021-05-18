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

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <vector>

namespace dhist {

template <typename T, bool kThreadsafe>
class InsertionBuffer;

template <typename T, bool kThreadsafe>
struct FlushIterator {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = value_type*;
  using reference = value_type&;

  FlushIterator(InsertionBuffer<T, kThreadsafe>* insertion_buffer)
      : ibuf_(insertion_buffer) {
    printf("> FlushIterator(): %zu, %zu\n", ibuf_->insert_buffer_begin_,
           ibuf_->insert_buffer_end_);

    ibuf_->flush_mu_.lock();

    idx_ = ibuf_->insert_buffer_begin_;
    ptr_ = &ibuf_->insert_buffer_[idx_];
  }

  FlushIterator(pointer ptr) : ibuf_(nullptr), ptr_(ptr) {}

  ~FlushIterator() {
    if (ibuf_) {
      ibuf_->insert_buffer_begin_ = idx_;

      ibuf_->flush_mu_.unlock();
    }
  }

  reference operator*() const { return *ptr_; }
  pointer operator->() { return ptr_; }

  // Prefix increment
  FlushIterator& operator++() {
    idx_++;
    if (idx_ >= ibuf_->insert_buffer_.size()) {
      idx_ = 0;
    }
    ptr_ = &ibuf_->insert_buffer_[idx_];
    return *this;
  }

  // Postfix increment
  FlushIterator operator++(int) {
    FlushIterator tmp = *this;
    ++(*this);
    return tmp;
  }

  friend bool operator==(const FlushIterator& a, const FlushIterator& b) {
    return a.ptr_ == b.ptr_;
  };

  friend bool operator!=(const FlushIterator& a, const FlushIterator& b) {
    return a.ptr_ != b.ptr_;
  };

 private:
  InsertionBuffer<T, kThreadsafe>* ibuf_;
  size_t idx_;
  pointer ptr_;
};

template <typename T, bool kThreadsafe = true>
class InsertionBuffer {
 public:
  friend struct FlushIterator<T, kThreadsafe>;

  InsertionBuffer(size_t buffer_size = 512)
      : generation_(0),
        refresh_generation_(0),
        insert_buffer_begin_(0),
        insert_buffer_end_(0) {
    insert_buffer_.resize(buffer_size);
  }

  void addValue(T val) {
    std::unique_ptr<std::scoped_lock<std::mutex>> lp;
    if (kThreadsafe) {
      lp.reset(new std::scoped_lock(flush_mu_));
    }

    insert_buffer_[insert_buffer_end_] = val;
    if (insert_buffer_end_ + 1 == insert_buffer_.size()) {
      insert_buffer_end_ = 0;
    } else {
      insert_buffer_end_++;
    }
  }

  FlushIterator<T, kThreadsafe> begin() {
    return FlushIterator(this);
  }

  FlushIterator<T, kThreadsafe> end() {
    return FlushIterator<T, kThreadsafe>(&insert_buffer_[insert_buffer_end_]);
  }

 protected:
  std::mutex flush_mu_;

 private:
  std::mutex insert_mu_;
  uint64_t generation_;
  uint64_t refresh_generation_;
  size_t insert_buffer_begin_;
  size_t insert_buffer_end_;
  std::vector<T> insert_buffer_;
};

}  // namespace dhist
