/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2023 Tobias Heuer <tobias.heuer@kit.edu>
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

#include "kahypar-resources/meta/policy_registry.h"

#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/macros.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/partition/context.h"

namespace mt_kahypar {

/**
 * The gain cache stores the gain values for all possible node moves for the cut metric.
 *
 * For a weighted hypergraph H = (V,E,c,w), the cut metric is defined as follows
 * connectivity(H) := \sum_{e \in cut(E)} w(e).
 *
 * The gain of moving a node u from its current block V_i to a target block V_j can be expressed as follows
 * g(u, V_j) := g(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) = |e| - 1 }) - w({ e \in I(u) | pin_count(e, V_i) = |e| }).
 * Moving node u from V_i to V_j, removes all nets e \in I(u) from the cut where pin_cout(e, V_j) = |e| - 1,
 * but makes it a cut hyperedge if pin_count(e, V_i) = |e|.
 *
 * We call the first term in the equation the benefit term b(u, V_j) and the second the penalty term p(u).
 * Our gain cache stores and maintains these entries for each node and block.
 * Thus, the gain cache stores k + 1 entries per node.
 * 
 * 
 * gain = benefit - penalty 
 * gain > 0 --> improvement
 * Dummy: benefit = -1, penalty = 1 --> never improvement
*/
class DummyGainCache {

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

  using AdjacentBlocksIterator = IntegerRangeIterator<PartitionID>::const_iterator;

 public:

  static constexpr GainPolicy TYPE = GainPolicy::conductance_local;
  static constexpr bool requires_notification_before_update = false;
  static constexpr bool initializes_gain_cache_entry_after_batch_uncontractions = false;
  static constexpr bool invalidates_entries = false;

  DummyGainCache() :
    _is_initialized(false),
    _k(kInvalidPartition),
    _dummy_adjacent_blocks() { }

  DummyGainCache(const Context&) :
    _is_initialized(false),
    _k(kInvalidPartition),
    _dummy_adjacent_blocks() { }

  DummyGainCache(const DummyGainCache&) = delete;
  DummyGainCache & operator= (const DummyGainCache &) = delete;

  DummyGainCache(DummyGainCache&& other) = default;
  DummyGainCache & operator= (DummyGainCache&& other) = default;

  // ####################### Initialization #######################

  bool isInitialized() const {
    return _is_initialized;
  }

  void reset(const bool run_parallel = true) {
    unused(run_parallel);
    _is_initialized = false;
  }

  size_t size() const {
    return 0;
  }

  // ! Initializes all gain cache entries
  template<typename PartitionedHypergraph>
  void initializeGainCache(const PartitionedHypergraph& partitioned_hg);

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void initializeGainCacheEntryForNode(const PartitionedHypergraph&,
                                       const HypernodeID&) {
    // Do nothing
  }

  IteratorRange<AdjacentBlocksIterator> adjacentBlocks(const HypernodeID) const {
    // We do not maintain the adjacent blocks of a node in this gain cache.
    // We therefore return an iterator over all blocks here
    return IteratorRange<AdjacentBlocksIterator>(
      _dummy_adjacent_blocks.cbegin(), _dummy_adjacent_blocks.cend());
  }

  // ####################### Gain Computation #######################

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, V_i) = |e| })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID /* only relevant for graphs */) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    unused(u);
    return 1;
  }

  // ! Recomputes the penalty term entry in the gain cache
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void recomputeInvalidTerms(const PartitionedHypergraph& partitioned_hg,
                             const HypernodeID u) {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    unused(partitioned_hg);
    unused(u);
    // Do nothing
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) = |e| - 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    unused(u);
    unused(to);
    return -1;
  }

  // ! Returns the gain of moving node u from its current block to a target block V_j.
  // ! More formally, g(u, V_j) := b(u, V_j) - p(u).
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID, /* only relevant for graphs */
                       const PartitionID to ) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    return benefitTerm(u, to) - penaltyTerm(u, kInvalidPartition);
  }

  // ####################### Delta Gain Update #######################

  // ! This function returns true if the corresponding syncronized edge update triggers
  // ! a gain cache update.
  static bool triggersDeltaGainUpdate(const SynchronizedEdgeUpdate& sync_update);

  // ! The partitioned (hyper)graph call this function when its updates its internal
  // ! data structures before calling the delta gain update function. The partitioned
  // ! (hyper)graph holds a lock for the corresponding (hyper)edge when calling this
  // ! function. Thus, it is guaranteed that no other thread will modify the hyperedge.
  template<typename PartitionedHypergraph>
  void notifyBeforeDeltaGainUpdate(const PartitionedHypergraph&, const SynchronizedEdgeUpdate&) {
    // Do nothing
  }

  // ! This functions implements the delta gain updates for the cut metric.
  // ! When moving a node from its current block from to a target block to, we iterate
  // ! over its incident hyperedges and update their pin count values. After each pin count
  // ! update, we call this function to update the gain cache to changes associated with
  // ! corresponding hyperedge.
  template<typename PartitionedHypergraph>
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const SynchronizedEdgeUpdate& sync_update);

  // ####################### Uncontraction #######################

  // ! This function implements the gain cache update after an uncontraction that restores node v in
  // ! hyperedge he. After the uncontraction node u and v are contained in hyperedge he.
  template<typename PartitionedHypergraph>
  void uncontractUpdateAfterRestore(const PartitionedHypergraph& partitioned_hg,
                                    const HypernodeID u,
                                    const HypernodeID v,
                                    const HyperedgeID he,
                                    const HypernodeID pin_count_in_part_after);

  // ! This function implements the gain cache update after an uncontraction that replaces u with v in
  // ! hyperedge he. After the uncontraction only node v is contained in hyperedge he.
  template<typename PartitionedHypergraph>
  void uncontractUpdateAfterReplacement(const PartitionedHypergraph& partitioned_hg,
                                        const HypernodeID u,
                                        const HypernodeID v,
                                        const HyperedgeID he);

  // ! This function is called after restoring a single-pin hyperedge. The function assumes that
  // ! u is the only pin of the corresponding hyperedge, while block_of_u is its corresponding block ID.
  void restoreSinglePinHyperedge(const HypernodeID,
                                 const PartitionID,
                                 const HyperedgeWeight) {
    // Do nothing here
  }

  // ! This function is called after restoring a net that became identical to another due to a contraction.
  template<typename PartitionedHypergraph>
  void restoreIdenticalHyperedge(const PartitionedHypergraph&,
                                 const HyperedgeID) {
    // Do nothing
  }

  // ! Notifies the gain cache that all uncontractions of the current batch are completed.
  void batchUncontractionsCompleted() {
    // Do nothing
  }

  // ####################### Only for Testing #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputePenaltyTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u) const {
    ASSERT(_is_initialized, "Gain cache is not initialized");
    unused(partitioned_hg);
    unused(u);
    HyperedgeWeight penalty = 1; 
    /* gain = benefit - penalty 
     * gain > 0 --> improvement
     * Dummy: benefit = -1, penalty = 1 --> never improvement
     */
    return penalty;
  }

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight recomputeBenefitTerm(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       const PartitionID to) const {
    unused(partitioned_hg);
    unused(u);
    unused(to);
    HyperedgeWeight benefit = -1; // benefit = old - new
    return benefit;
  }

  void changeNumberOfBlocks(const PartitionID new_k) {
    ASSERT(new_k <= _k);
    _dummy_adjacent_blocks = IntegerRangeIterator<PartitionID>(new_k);
  }

  template<typename PartitionedHypergraph>
  bool verifyTrackedAdjacentBlocksOfNodes(const PartitionedHypergraph&) const {
    // Gain cache does not track adjacent blocks of node
    return true;
  }

 private:
  friend class DeltaDummyGainCache;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t penalty_index(const HypernodeID u) const {
    unused(u);
    return 0;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  size_t benefit_index(const HypernodeID u, const PartitionID p) const {
    unused(u);
    unused(p);
    return 0;
  }

  // ! Allocates the memory required to store the gain cache
  void allocateGainTable(const HypernodeID num_nodes,
                         const PartitionID k) {
    unused(num_nodes);
    if (k != kInvalidPartition) {
      _k = k;
      _dummy_adjacent_blocks = IntegerRangeIterator<PartitionID>(k);
    }
  }

  // ! Initializes the benefit and penalty terms for a node u
  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                       const HypernodeID u,
                                       vec<Gain>& benefit_aggregator);

  bool nodeGainAssertions(const HypernodeID u, const PartitionID p) const {
    if ( p == kInvalidPartition || p >= _k ) {
      LOG << "Invalid block ID (Node" << u << "is part of block" << p
          << ", but valid block IDs must be in the range [ 0," << _k << "])";
      return false;
    }
    return true;
  }


  // ! Indicate whether or not the gain cache is initialized
  bool _is_initialized;

  // ! Number of blocks
  PartitionID _k;

  // ! Provides an iterator from 0 to k (:= number of blocks)
  IntegerRangeIterator<PartitionID> _dummy_adjacent_blocks;
};

/**
 * In our FM algorithm, the different local searches perform nodes moves locally not visible for other
 * threads. The delta gain cache stores these local changes relative to the shared
 * gain cache. For example, the penalty term can be computed as follows
 * p'(u) := p(u) + Δp(u)
 * where p(u) is the penalty term stored in the shared gain cache and Δp(u) is the penalty term stored in
 * the delta gain cache after performing some moves locally. To maintain Δp(u) and Δb(u,V_j), we use a hash
 * table that only stores entries affected by a gain cache update.
*/
class DeltaDummyGainCache {

  using AdjacentBlocksIterator = typename DummyGainCache::AdjacentBlocksIterator;

 public:
  static constexpr bool requires_connectivity_set = false;

  DeltaDummyGainCache(const DummyGainCache& gain_cache) :
    _gain_cache(gain_cache)  { }

  // ####################### Initialize & Reset #######################

  void initialize(const size_t size) {
    unused(size);
  }

  void clear() {
    // do nothing
  }

  void dropMemory() {
    // do nothing
  }

  size_t size_in_bytes() const {
    return 0;
  }

  // ####################### Gain Computation #######################

  // ! Returns an iterator over the adjacent blocks of a node
  IteratorRange<AdjacentBlocksIterator> adjacentBlocks(const HypernodeID hn) const {
    return _gain_cache.adjacentBlocks(hn);
  }

  // ! Returns the penalty term of node u.
  // ! More formally, p(u) := (w(I(u)) - w({ e \in I(u) | pin_count(e, V_i) = |e| })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight penaltyTerm(const HypernodeID u,
                              const PartitionID from) const {
    return _gain_cache.penaltyTerm(u, from);
  }

  // ! Returns the benefit term for moving node u to block to.
  // ! More formally, b(u, V_j) := w({ e \in I(u) | pin_count(e, V_j) = |e| - 1 })
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight benefitTerm(const HypernodeID u, const PartitionID to) const {
    ASSERT(to != kInvalidPartition && to < _gain_cache._k);
    return _gain_cache.benefitTerm(u, to);
  }

  // ! Returns the gain of moving node u from its current block to a target block V_j.
  // ! More formally, g(u, V_j) := b(u, V_j) - p(u).
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight gain(const HypernodeID u,
                       const PartitionID from,
                       const PartitionID to ) const {
    return benefitTerm(u, to) - penaltyTerm(u, from);
  }

 // ####################### Delta Gain Update #######################

  template<typename PartitionedHypergraph>
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  void deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                       const SynchronizedEdgeUpdate& sync_update) {
    unused(partitioned_hg);
    unused(sync_update);
    // do nothing
  }

 // ####################### Miscellaneous #######################

  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    unused(parent);
    // utils::MemoryTreeNode* gain_cache_delta_node = parent->addChild("Delta Gain Cache");
    // gain_cache_delta_node->updateSize(size_in_bytes());
  }

 private:
  const DummyGainCache& _gain_cache;
};

}  // namespace mt_kahypar
