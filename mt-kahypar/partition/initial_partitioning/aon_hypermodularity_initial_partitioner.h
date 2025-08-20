/*******************************************************************************
 * MIT License
 * 
 * This file is part of Mt-KaHyPar.
 *
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

#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/array.h"

namespace mt_kahypar {

template<typename TypeTraits>
class AONHypermodularityInitialPartitioner : public IInitialPartitioner {

  static constexpr bool debug = false;

  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using UnderlyingHypergraph = typename PartitionedHypergraph::UnderlyingHypergraph;

 public:
  AONHypermodularityInitialPartitioner(const InitialPartitioningAlgorithm,
                           ip_data_container_t* ip_data,
                           const Context& context,
                           const int seed, const int tag) :
    _ip_data(ip::to_reference<TypeTraits>(ip_data)),
    _context(context),
    _rng(seed),
    _tag(tag) { }

 private:
  void partitionImpl() final;

  bool fitsIntoBlock(PartitionedHypergraph& hypergraph,
                     const HypernodeID hn,
                     const PartitionID block) const {
    ASSERT(block != kInvalidPartition && block < _context.partition.k);
    return hypergraph.partWeight(block) + hypergraph.nodeWeight(hn) <=
      _context.partition.perfect_balance_part_weights[block];
  }

  // Used if some vertices have fixed labels
  inline void randomPartition(PartitionedHypergraph& hypergraph);

  // Contract communities of the coarsest hypergraph and rewrites its partition
  inline void collapse(UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z);

  // Perform the Louvain step on the collapsed hypergraph
  void louvainStep(UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, const vec<double>& beta, const vec<double>& gamma, const long long maxNumIter=100, const double eps=1e-8, const bool randomize=true);

  // Perform the Louvain step for a single node on the collapsed hypergraph
  inline bool louvainStepForANode(const HypernodeID& i, const vec<HypernodeID>& neighbors_i, const ds::Array<bool>& visitedParts, UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, const vec<double>& beta, const vec<double>& gamma, const long long maxNumIter=100, const double eps=1e-8, const bool randomize=true);

  // Calculate the gain of moving node i to partition A
  // using the AllOrNothing-Hypermodularity gain function
  inline double QAONGain(PartitionedHypergraph& H_new_partitioned, const HypernodeID i, const PartitionID A, const vec<double>& beta, const vec<double>& gamma);

  // Adjust current communities and check if they changed
  inline bool expand(UnderlyingHypergraph& H, UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, vec<HypernodeID>& z);

  InitialPartitioningDataContainer<TypeTraits>& _ip_data;
  const Context& _context;
  std::mt19937 _rng;
  const int _tag;
};

} // namespace mt_kahypar
