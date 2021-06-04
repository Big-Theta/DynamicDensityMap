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

#include <atomic>
#include <cstddef>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

namespace dyden {

template <typename T>
class InsertionBuffer;

template <typename T>
struct FlushIterator {
  using iterator_category = std::forward_iterator_tag;
  using difference_type = std::ptrdiff_t;
  using value_type = T;
  using pointer = value_type*;
  using reference = value_type&;

  FlushIterator(InsertionBuffer<T>* insertion_buffer)
      : ibuf_(insertion_buffer), processed_(0) {
    ibuf_->flush_mu_.lock();
    idx_ = ibuf_->buffer_begin_;
    end_idx_ = ibuf_->buffer_end_;
    ptr_ = &ibuf_->buffer_[idx_];
  }

  FlushIterator(pointer ptr) : ibuf_(nullptr), ptr_(ptr) {}

  FlushIterator(const FlushIterator&) = delete;

  ~FlushIterator() {
    if (ibuf_) {
      ibuf_->buffer_begin_ = idx_;

      std::scoped_lock(ibuf_->insert_mu_);
      ibuf_->unflushed_ -= processed_;
      ibuf_->flush_mu_.unlock();
    }
  }

  reference operator*() const { return *ptr_; }
  pointer operator->() { return ptr_; }

  operator bool() const {
    return idx_ != end_idx_;
  }

  // Prefix increment
  FlushIterator& operator++() {
    processed_++;
    idx_++;
    if (idx_ >= ibuf_->buffer_.size()) {
      idx_ = 0;
    }
    ptr_ = &ibuf_->buffer_[idx_];
    return *this;
  }

  friend bool operator==(const FlushIterator& a, const FlushIterator& b) {
    return a.ptr_ == b.ptr_;
  };

  friend bool operator!=(const FlushIterator& a, const FlushIterator& b) {
    return a.ptr_ != b.ptr_;
  };

 private:
  InsertionBuffer<T>* ibuf_;
  size_t idx_;
  size_t end_idx_;
  size_t processed_;
  pointer ptr_;
};

template <typename T>
class InsertionBuffer {
 public:
  friend struct FlushIterator<T>;

  InsertionBuffer(size_t buffer_size = 512)
      : buffer_begin_(0), buffer_end_(0), unflushed_(0) {
    buffer_.resize(buffer_size);
  }

  size_t addValue(T val) {
    std::scoped_lock l(insert_mu_);

    buffer_[buffer_end_] = val;
    size_t end = 1 + buffer_end_;
    if (end == buffer_.size()) {
      end = 0;
    }
    buffer_end_ = end;

    unflushed_++;
    return unflushed_;
  }

  size_t capacity() const {
    return buffer_.size();
  }

  FlushIterator<T> lockedIterator() {
    return FlushIterator(this);
  }

 protected:
  std::mutex insert_mu_;
  std::mutex flush_mu_;

 private:
  size_t buffer_begin_;
  size_t buffer_end_;
  size_t unflushed_;
  std::vector<T> buffer_;
};

}  // namespace dyden
