// MIT License
//
// Copyright (c) 2021 Logan Evans
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

#include <assert.h>

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

namespace dyden {

template <typename T>
class LocklessInsertionBuffer;

template <typename T>
struct LockedFlushIterator {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = value_type*;
  using reference = value_type&;

  LockedFlushIterator(LocklessInsertionBuffer<T>* insertion_buffer)
      : ibuf_(insertion_buffer), index_mask_(ibuf_->index_mask_) {
    ibuf_->flush_mu_.lock();

    num_inserted_ = ibuf_->num_inserted_.load(std::memory_order_acquire);
    // Writers must have flush_mu_, so there is no data race with reading
    // num_removed_ here.
    num_removed_ = ibuf_->num_removed_.load(std::memory_order_relaxed);

    // Check if there are dropped entries.
    if (num_inserted_ - num_removed_ > ibuf_->capacity()) {
      num_removed_ = num_inserted_ - ibuf_->capacity();
    }

    // Permit the use of ibuf_->buffer_[X].load(std::memory_order_relaxed).
    std::atomic_thread_fence(std::memory_order_acquire);

    //current_val_ =
    //    ibuf_->buffer_[index(num_removed_)].load(std::memory_order_relaxed);
    //num_removed_++;

    if (num_removed_ < num_inserted_) {
      size_t idx = index(num_removed_);
      current_val_ = ibuf_->buffer_[idx].load(std::memory_order_relaxed);
      // TODO: Use a template option for this sentinel so that it works with
      // std::pair<double, double>.
      ibuf_->buffer_[idx].store(std::numeric_limits<T>::max());
    }
    num_removed_++;
  }

  LockedFlushIterator(pointer ptr) : ibuf_(nullptr), num_removed_(0) {}

  LockedFlushIterator(const LockedFlushIterator&) = delete;

  ~LockedFlushIterator() {
    if (ibuf_) {
      // Flush each ibuf_->buffer_[X].store(V, std::memory_order_relaxed).
      std::atomic_thread_fence(std::memory_order_release);
      ibuf_->num_removed_.store(num_inserted_, std::memory_order_release);
      ibuf_->flush_mu_.unlock();
    }
  }

  value_type operator*() const {
    return current_val_;
  }
  pointer operator->() { return &current_val_; }

  operator bool() const {
    bool val = num_removed_ <= num_inserted_;
    return val;
  }

  // Prefix increment
  LockedFlushIterator& operator++() {
    size_t idx;

    do {
      idx = index(num_removed_);
      num_removed_++;
      current_val_ = ibuf_->buffer_[idx].load(std::memory_order_relaxed);
      // This check guards against a situation where a separate thread has
      // claimed an index in the buffer but has not yet inserted its value.
      // The value will be dropped.
    } while (num_removed_ <= num_inserted_ &&
             current_val_ == std::numeric_limits<T>::max());

    if (num_removed_ <= num_inserted_) {
      // TODO: Use a template option for this sentinel so that it works with
      // std::pair<double, double>.
      ibuf_->buffer_[idx].store(std::numeric_limits<T>::max());
    }

    return *this;
  }

  friend bool operator==(const LockedFlushIterator& a,
                         const LockedFlushIterator& b) {
    return a.ibuf_ == b.ibuf_ && a.num_removed_ == b.num_removed_;
  };

  friend bool operator!=(const LockedFlushIterator& a,
                         const LockedFlushIterator& b) {
    return a.num_removed_ != b.num_removed_ || a.ibuf_ != b.ibuf_;
  };

 private:
  LocklessInsertionBuffer<T>* ibuf_;
  const size_t index_mask_;
  size_t num_inserted_;
  size_t num_removed_;
  value_type current_val_;

  size_t index(size_t num_removed) const { return num_removed & index_mask_; }
};


template <typename T>
class LocklessInsertionBuffer {
 public:
  friend struct LockedFlushIterator<T>;

  LocklessInsertionBuffer(size_t buffer_size = 512)
      : buffer_size_(buffer_size),
        index_mask_(buffer_size_ - 1),
        flush_at_population_(buffer_size_ / 2),
        num_inserted_(0),
        num_removed_(0) {
    assert((buffer_size & index_mask_) == 0);
    buffer_ = new std::atomic<T>[buffer_size_];
  }

  ~LocklessInsertionBuffer() {
    delete[] buffer_;
  }

  // Returns bool: true if the caller should flush the queue.
  bool addValue(T val) {
    size_t inserted = num_inserted_.fetch_add(1, std::memory_order_acq_rel);
    size_t idx = index(inserted);
    buffer_[idx].store(val, std::memory_order_release);

    size_t population = inserted - num_removed_.load(std::memory_order_acquire);
    if (population < flush_at_population_) {
      return false;
    }

    if (population == flush_at_population_) {
      return true;
    }

    if (population + 2 < buffer_size_) {
      return false;
    }

    // Apply some backpressure so that we can reduce the amount of dropped
    // entries.
    std::scoped_lock l(flush_mu_);
    return true;
  }

  size_t capacity() const { return buffer_size_; }

  LockedFlushIterator<T> lockedIterator() {
    return LockedFlushIterator(this);
  }

 private:
  const size_t buffer_size_;
  const size_t index_mask_;
  const size_t flush_at_population_;
  std::mutex flush_mu_;
  std::atomic<size_t> num_inserted_;
  std::atomic<size_t> num_removed_;
  std::atomic<T>* buffer_;

  size_t index(size_t num_inserted) const { return num_inserted & index_mask_; }
};

}  // namespace dyden
