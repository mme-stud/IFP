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

#include <vector>
 
#include <tbb/enumerable_thread_specific.h>


#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include "mt-kahypar/partition/refinement/gains/conductance_global/conductance_global_attributed_gains.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/conductance_pq.h"
#include "mt-kahypar/datastructures/static_bitset.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
 
class ConductanceLocalGainComputation : public GainComputationBase<ConductanceLocalGainComputation, ConductanceGlobalAttributedGains> {
  using Base = GainComputationBase<ConductanceLocalGainComputation, ConductanceGlobalAttributedGains>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;
  static constexpr size_t BITS_PER_BLOCK = ds::StaticBitset::BITS_PER_BLOCK;

public:
ConductanceLocalGainComputation(const Context& context,
                               bool disable_randomization = false) :
    // [took from steiner_tree_gain_computation]
    Base(context, disable_randomization),
    _local_adjacent_blocks([&] { return constructBitset(); }),
    _all_blocks(context.partition.k) {
    for ( PartitionID to = 0; to < context.partition.k; ++to )  {
      _all_blocks.set(to);
    }
  }
 
  // ! Precomputes the gain to all adjacent blocks.
  // ! Conceptually, we compute the gain of moving the node to an non-adjacent block
  // ! and the gain to all adjacent blocks assuming the node is in an isolated block.
  // ! The gain of that node to a block to can then be computed by
  // ! 'isolated_block_gain - tmp_scores[to]' (see gain(...))
  // ! (new) Correction: 
  // ! The gain to non-adjacent blocks could be different for conductance
  // ! => put in tmp_scores the gain to considered blocks
  // !    calculate the gain = tmp_scores[to]
  // !    set isolated_block_gain = 0
  template<typename PartitionedHypergraph>
  void precomputeGains(const PartitionedHypergraph& phg,
                       const HypernodeID hn,
                       RatingMap& tmp_scores,
                       Gain& isolated_block_gain,
                       const bool consider_non_adjacent_blocks) {
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");
    isolated_block_gain = 0;

    // Compute all adjacent blocks of node 
    // [took from steiner_tree_gain_computation.h]
    ds::Bitset& adjacent_blocks = consider_non_adjacent_blocks ?
                                  _all_blocks : _local_adjacent_blocks.local();
    ds::StaticBitset adjacent_blocks_view(
                          adjacent_blocks.numBlocks(), adjacent_blocks.data());
    if ( !consider_non_adjacent_blocks ) {
      adjacent_blocks.reset();
      for (const HyperedgeID& he : phg.incidentEdges(hn)) {
        for ( const PartitionID& to : phg.connectivitySet(he) ) {
          adjacent_blocks.set(to);
        }
      }
      adjacent_blocks.unset(phg.partID(hn));
    }
    
    // Set the trivial values needed for emulating move
    PartitionID from = phg.partID(hn);
    HypergraphVolume cut_weight_from_before = phg.partCutWeight(from);
    bool original_stats = phg.conductancePriorityQueueUsesOriginalStats();
    HypergraphVolume volume_from_before = 0;
    HypergraphVolume volume_from_after = 0;
    HypergraphVolume total_volume = 0;
    HypergraphVolume weighted_degree = 0;
    PartitionID k = _context.partition.k;
    if (original_stats) {
      weighted_degree = phg.nodeOriginalWeightedDegree(hn);
      total_volume = phg.originalTotalVolume();
      volume_from_before = phg.partOriginalVolume(from);
    } else {
      weighted_degree = phg.nodeWeightedDegree(hn);
      total_volume = phg.totalVolume();
      volume_from_before = phg.partVolume(from);
    }
    volume_from_after = volume_from_before - weighted_degree;

    // Current conductance objective of the "from" part
    ASSERT(volume_from_before <= total_volume, "Volume of partition " << V(from) << " is larger than total volume: " << V(volume_from_before) << " > " << V(total_volume));
    ds::ConductanceFraction fraction_from_before(cut_weight_from_before, std::min(volume_from_before, total_volume - volume_from_before));

    // Compute the gain to all concident blocks
    for ( const PartitionID to : adjacent_blocks_view )  {
      HypergraphVolume volume_to_before = 0;
      HypergraphVolume volume_to_after = 0;
      HypergraphVolume cut_weight_from_after = 0;
      HypergraphVolume cut_weight_to_before = 0;
      HypergraphVolume cut_weight_to_after = 0;
      // sets the right versions of cut weights as before the move, 
      // the right version of volume_to as after the move
      setPreIncidentEdgeLoopData(phg, from, to, original_stats, weighted_degree,
                                 cut_weight_from_before, cut_weight_to_before,
                                 cut_weight_from_after, cut_weight_to_after,
                                 volume_to_before, volume_to_after);

      // Current conductance objective of the "to" part
      ASSERT(volume_to_before < total_volume, "Volume of partition " << V(to) << " is larger or equal to total volume: " << V(volume_to_before) << " > " << V(total_volume));
      ds::ConductanceFraction fraction_to_before(cut_weight_to_before, std::min(volume_to_before, total_volume - volume_to_before));
      
      // Current LOCAL conductance objective between "from" and "to" parts
      ds::ConductanceFraction max_fraction_before = fraction_to_before;
      if (fraction_from_before > fraction_to_before) {
        max_fraction_before = fraction_from_before;
      }
      HyperedgeWeight local_conductance_before = 
            ConductanceGlobalAttributedGains::compute_conductance_objective(
              total_volume, max_fraction_before, k); // at least one part wasn't empty

      // adjust cut weights
      for (HyperedgeID he : phg.incidentEdges(hn)) {
        HyperedgeID edge_size = phg.edgeSize(he);
        if (edge_size == 1) {
          // Never a cutting edge
          continue;
        }
        HypernodeID pin_count_in_from_part = phg.pinCountInPart(he, from);
        HypernodeID pin_count_in_to_part = phg.pinCountInPart(he, to);
        HyperedgeWeight he_weight = phg.edgeWeight(he);
        if (pin_count_in_from_part == 1) {
          // Was a cutting edge before for "from" part
          // but not anymore
          cut_weight_from_after -= he_weight;
        } else if (pin_count_in_from_part == edge_size) {
          // Wasn't a cutting edge before for "from" part
          // but is now
          cut_weight_from_after += he_weight;
        }
        if (pin_count_in_to_part == 0) {
          // Wasn't a cutting edge before for "to" part
          // but is now
          cut_weight_to_after += he_weight;
        } else if (pin_count_in_to_part == edge_size - 1) {
          // Was a cutting edge before for "to" part
          // but not anymore
          cut_weight_to_after -= he_weight;
        }
      }

      // Get the new conductance objective of the "from" part
      ASSERT(volume_from_after < total_volume, "Volume of partition " << V(from) << " is larger or equal to total volume: " << V(volume_from_after) << " > " << V(total_volume));
      ds::ConductanceFraction fraction_from_after(cut_weight_from_after, std::min(volume_from_after, total_volume - volume_from_after));
      // Get the new conductance objective of the "to" part
      ASSERT(volume_to_after <= total_volume, "Volume of partition " << V(to) << " is larger than total volume: " << V(volume_to_after) << " > " << V(total_volume));
      ds::ConductanceFraction fraction_to_after(cut_weight_to_after, std::min(volume_to_after, total_volume - volume_to_after));

      // Get the new conductance objective      
      ds::ConductanceFraction max_fraction_after = fraction_to_after;
      if (fraction_from_after > fraction_to_after) {
        max_fraction_after = fraction_from_after;
      }
      HyperedgeWeight local_conductance_after = 
            ConductanceGlobalAttributedGains::compute_conductance_objective(
              total_volume, max_fraction_after, k); // at least one part isn't empty

      tmp_scores[to] = local_conductance_after - local_conductance_before;
    }
  }
 
  HyperedgeWeight gain(const Gain to_score,
                       const Gain) {
    return to_score;
  }

  void changeNumberOfBlocksImpl(const PartitionID new_k) {
    ASSERT(new_k == _context.partition.k);
    // [took from steiner_tree_gain_computation]
    for ( auto& adjacent_blocks : _local_adjacent_blocks ) {
      adjacent_blocks.resize(new_k);
    }
    _all_blocks.resize(new_k);
    for ( PartitionID to = 0; to < new_k; ++to )  {
      _all_blocks.set(to);
    }
  }

 private:
  ds::Bitset constructBitset() const {
    return ds::Bitset(_context.partition.k);
  }

  // ! This function is used to set almost all the data needed
  // ! for the emulation of move:
  // ! - cut_weight_from_before: the cut weight of the partition from before the move
  // ! - cut_weight_to_before: the cut weight of the partition to before the move
  // ! - cut_weight_from_after: the cut weight of the partition from after the move
  // ! - cut_weight_to_after: the cut weight of the partition to after the move
  // ! - volume_to_before: the used version of volume of the partition to before the move
  // ! - volume_to_after: the used version of volume of the partition to after the move
  // ! Note: cut weights have to be adjusted later in the precomputeGains function
  template<typename PartitionedHypergraph>
  void setPreIncidentEdgeLoopData(const PartitionedHypergraph& phg,
                                const PartitionID from,
                                const PartitionID to,
                                const bool original_stats,
                                const HypergraphVolume weighted_degree,
                                HypergraphVolume& cut_weight_from_before,
                                HypergraphVolume& cut_weight_to_before,
                                HypergraphVolume& cut_weight_from_after,
                                HypergraphVolume& cut_weight_to_after,
                                HypergraphVolume& volume_to_before,
                                HypergraphVolume& volume_to_after) {
    ASSERT(0 <= to && to <= _context.partition.k, "Invalid partition ID: " << V(to) << " k=" << _context.partition.k);
    ASSERT(0 <= from && from <= _context.partition.k, "Invalid partition ID: " << V(from) << " k=" << _context.partition.k);
    ASSERT(from != to, "Invalid partition ID: from=to=" << from);
    cut_weight_from_before = phg.partCutWeight(from);
    cut_weight_to_before = phg.partCutWeight(to);
    cut_weight_from_after = cut_weight_from_before; // will be adjusted
    cut_weight_to_after = cut_weight_to_before; // will be adjusted
    if (original_stats) {
      volume_to_before = phg.partOriginalVolume(to);
    } else {
      volume_to_before = phg.partVolume(to);
    }
    volume_to_after = volume_to_before + weighted_degree;
  }

  using Base::_context;

  // ! Before gain computation, we construct a bitset that contains all
  // ! adjacent nodes of a block
  // [took from steiner_tree_gain_computation]
  tbb::enumerable_thread_specific<ds::Bitset> _local_adjacent_blocks;
  ds::Bitset _all_blocks;
};
 
}  // namespace mt_kahypar