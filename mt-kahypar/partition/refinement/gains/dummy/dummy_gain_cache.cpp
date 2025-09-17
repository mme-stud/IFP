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

#include "mt-kahypar/partition/refinement/gains/dummy/dummy_gain_cache.h"

#include <tbb/parallel_for.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/concurrent_vector.h>

#include "mt-kahypar/definitions.h"

namespace mt_kahypar {

template<typename PartitionedHypergraph>
void DummyGainCache::initializeGainCache(const PartitionedHypergraph& partitioned_hg) {
  ASSERT(!_is_initialized, "Gain cache is already initialized");
  ASSERT(_k <= 0 || _k >= partitioned_hg.k(),
    "Gain cache was already initialized for a different k" << V(_k) << V(partitioned_hg.k()));
  allocateGainTable(partitioned_hg.topLevelNumNodes(), partitioned_hg.k());
  _is_initialized = true;
}

bool DummyGainCache::triggersDeltaGainUpdate(const SynchronizedEdgeUpdate& sync_update) {
  unused(sync_update);
  return false;
}


template<typename PartitionedHypergraph>
void DummyGainCache::deltaGainUpdate(const PartitionedHypergraph& partitioned_hg,
                                   const SynchronizedEdgeUpdate& sync_update) {
  ASSERT(_is_initialized, "Gain cache is not initialized");
  // do nothing
  unused(partitioned_hg);
  unused(sync_update);
}

template<typename PartitionedHypergraph>
void DummyGainCache::uncontractUpdateAfterRestore(const PartitionedHypergraph& partitioned_hg,
                                                const HypernodeID u,
                                                const HypernodeID v,
                                                const HyperedgeID he,
                                                const HypernodeID pin_count_in_part_after) {
  unused(partitioned_hg);
  unused(u);
  unused(v);
  unused(he);
  unused(pin_count_in_part_after);
  // do nothing
}

template<typename PartitionedHypergraph>
void DummyGainCache::uncontractUpdateAfterReplacement(const PartitionedHypergraph& partitioned_hg,
                                                    const HypernodeID u,
                                                    const HypernodeID v,
                                                    const HyperedgeID he) {
  unused(partitioned_hg);
  unused(u);
  unused(v);
  unused(he);
  // do nothing
}

template<typename PartitionedHypergraph>
void DummyGainCache::initializeGainCacheEntryForNode(const PartitionedHypergraph& partitioned_hg,
                                                   const HypernodeID u,
                                                   vec<Gain>& benefit_aggregator) {
  unused(partitioned_hg);
  unused(u);
  unused(benefit_aggregator);
  // do nothing
}


namespace {
#define DUMMY_INITIALIZE_GAIN_CACHE(X) void DummyGainCache::initializeGainCache(const X&)
#define DUMMY_DELTA_GAIN_UPDATE(X) void DummyGainCache::deltaGainUpdate(const X&,                     \
                                                                    const SynchronizedEdgeUpdate&)
#define DUMMY_RESTORE_UPDATE(X) void DummyGainCache::uncontractUpdateAfterRestore(const X&,          \
                                                                              const HypernodeID, \
                                                                              const HypernodeID, \
                                                                              const HyperedgeID, \
                                                                              const HypernodeID)
#define DUMMY_REPLACEMENT_UPDATE(X) void DummyGainCache::uncontractUpdateAfterReplacement(const X&,            \
                                                                                      const HypernodeID,   \
                                                                                      const HypernodeID,   \
                                                                                      const HyperedgeID)
#define DUMMY_INIT_GAIN_CACHE_ENTRY(X) void DummyGainCache::initializeGainCacheEntryForNode(const X&,           \
                                                                                        const HypernodeID,  \
                                                                                        vec<Gain>&)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DUMMY_INITIALIZE_GAIN_CACHE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DUMMY_DELTA_GAIN_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DUMMY_RESTORE_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DUMMY_REPLACEMENT_UPDATE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(DUMMY_INIT_GAIN_CACHE_ENTRY)

}  // namespace mt_kahypar
