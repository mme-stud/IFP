/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2019 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2019 Tobias Heuer <tobias.heuer@kit.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <cmath>
#include <vector>

#include "kahypar-resources/meta/mandatory.h"

#include "mt-kahypar/macros.h"


namespace mt_kahypar {
namespace ds {

/*!
 * Non concurrent version of ConcurrentBucketMap:
 * Non concurrent data structures to distribute key-value pairs
 * into buckets such that values with the same key reside in
 * the same bucket. Main use case for that data structure is
 * within the parallel hyperedge detection inside the AON_HMLL IP
 * contractions. There, hyperedges are inserted with its footprint
 * as key into this data structure and afterwards all hyperedges
 * with the same footprint reside in the same bucket. Finally,
 * parallel hyperedges can be detected by processing each bucket
 * sequential.
 * If we insert a key-value pair, we compute its corresponding
 * bucket by computing key % num_buckets. Note, key must be of 
 * type uint64_t.
 */
template <typename Value>
class SequentialBucketMap {

  static constexpr bool debug = false;
  static constexpr size_t BUCKET_FACTOR = 128;

  using Bucket = std::vector<Value>;

 public:

  SequentialBucketMap() :
    _num_buckets(align_to_next_power_of_two(BUCKET_FACTOR * 1)), // * std::thread::hardware_concurrency())),
    _mod_mask(_num_buckets - 1),
    _buckets(_num_buckets) { }

  SequentialBucketMap(const SequentialBucketMap&) = delete;
  SequentialBucketMap & operator= (const SequentialBucketMap &) = delete;

  SequentialBucketMap(SequentialBucketMap&& other) :
    _num_buckets(other._num_buckets),
    _mod_mask(_num_buckets - 1),
    _buckets(std::move(other._buckets)) { }

  // ! Applies function f to all buckets sequentially (!)
  // ! The name "doParallelForAllBuckets" is chosen to have the same
  // ! interface as ConcurrentBucketMap
  template<typename F>
  void doParallelForAllBuckets(const F& f) {
    for (size_t i = UL(0); i < _num_buckets; ++i) {
      f(i);
    }
  }

  // ! Returns the number of buckets
  size_t numBuckets() const {
    return _num_buckets;
  }

  // ! Returns the corresponding bucket
  Bucket& getBucket(const size_t bucket) {
    ASSERT(bucket < _num_buckets);
    return _buckets[bucket];
  }

  // ! Reserves memory in each bucket such that the estimated number of insertions
  // ! can be handled without the need (with high probability) of expensive bucket resizing.
  void reserve_for_estimated_number_of_insertions(const size_t estimated_num_insertions) {
    // ! Assumption is that keys are evenly distributed among buckets (with a small buffer)
    const size_t estimated_bucket_size = std::max(
      static_cast<size_t>( 1.5 * estimated_num_insertions ) / _num_buckets, UL(1));
    for (size_t i = 0; i < _num_buckets; ++i) {
      _buckets[i].reserve(estimated_bucket_size);
    }
  }

  // ! Inserts a key-value pair
  void insert(const size_t& key, Value&& value) {
    size_t bucket = key & _mod_mask;
    ASSERT(bucket < _num_buckets);
    _buckets[bucket].emplace_back( std::move(value) );
  }

  // ! Frees the memory of all buckets
  void free() {
     _buckets = std::vector<Bucket>();
  }

  // ! Frees the memory of the corresponding bucket
  void free(const size_t bucket) {
    ASSERT(bucket < _num_buckets);
    _buckets[bucket] = Bucket();
  }

  // ! Clears the corresponding bucket
  void clear(const size_t bucket) {
    ASSERT(bucket < _num_buckets);
    _buckets[bucket].clear();
  }

  // ! Clears all buckets sequentially (!)
  // ! The name "clearParallel" is chosen to have the same
  // ! interface as ConcurrentBucketMap
  void clearParallel() {
    doParallelForAllBuckets([&](const size_t i) {
      clear(i);
    });
  }

 private:
  size_t align_to_next_power_of_two(const size_t size) const {
    return std::pow(2.0, std::ceil(std::log2(static_cast<double>(size))));
  }

  const size_t _num_buckets;
  const size_t _mod_mask;
  std::vector<Bucket> _buckets;
};
}  // namespace ds
}  // namespace mt_kahypar
