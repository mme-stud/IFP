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
#include "mt-kahypar/datastructures/static_bitset.h"
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"

namespace mt_kahypar {
 
class ConductanceGlobalGainComputation : public GainComputationBase<ConductanceGlobalGainComputation, ConductanceGlobalAttributedGains> {
  using Base = GainComputationBase<ConductanceGlobalGainComputation, ConductanceGlobalAttributedGains>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;
  static constexpr size_t BITS_PER_BLOCK = ds::StaticBitset::BITS_PER_BLOCK;

public:
ConductanceGlobalGainComputation(const Context& context,
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
    }
    
    // Set the trivial values needed for emulating sync_update of conductance
    bool original_stats = phg.conductancePriorityQueueUsesOriginalStats();
    SynchronizedEdgeUpdate sync_update;
    sync_update.k = _context.partition.k;
    sync_update.from = phg.partID(hn);
    if (original_stats) {
      sync_update.weighted_degree = phg.nodeOriginalWeightedDegree(hn);
      sync_update.total_volume = phg.originalTotalVolume();
      sync_update.volume_from_after = phg.partOriginalVolume(sync_update.from) - sync_update.weighted_degree;
    } else {
      sync_update.weighted_degree = phg.nodeWeightedDegree(hn);
      sync_update.total_volume = phg.totalVolume();
      sync_update.volume_from_after = phg.partVolume(sync_update.from) - sync_update.weighted_degree;
    }
    sync_update.top_three_conductance_info_before = phg.topThreePartConductanceInfos();

    // Current conductance objective
    ds::ConductanceFraction conductance_fraction_now = phg.topPartConductanceInfo().fraction;
    HyperedgeWeight conductance_now = ConductanceGlobalAttributedGains::
          compute_conductance_objective(sync_update.total_volume, conductance_fraction_now, sync_update.k);

    // Compute the gain to all concident blocks
    for ( const PartitionID to : adjacent_blocks_view )  {
      // sets the right versions of cut weights as before the move, 
      // the right version of volume_to as after the move
      setPreSyncUpdataData(phg, to, original_stats, sync_update);

      // adjust cut weights
      for (HyperedgeID he : phg.incidentEdges(hn)) {
        HyperedgeID edge_size = phg.edgeSize(he);
        if (edge_size == 1) {
          // Never a cutting edge
          continue;
        }
        HypernodeID pin_count_in_from_part = phg.pinCountInPart(he, sync_update.from);
        HypernodeID pin_count_in_to_part = phg.pinCountInPart(he, to);
        HyperedgeWeight he_weight = phg.edgeWeight(he);
        if (pin_count_in_from_part == 1) {
          // Was a cutting edge before for "from" part
          // but not anymore
          sync_update.cut_weight_from_after -= he_weight;
        } else if (pin_count_in_from_part == edge_size) {
          // Wasn't a cutting edge before for "from" part
          // but is now
          sync_update.cut_weight_from_after += he_weight;
        }
        if (pin_count_in_to_part == 0) {
          // Wasn't a cutting edge before for "to" part
          // but is now
          sync_update.cut_weight_to_after += he_weight;
        } else if (pin_count_in_to_part == edge_size - 1) {
          // Was a cutting edge before for "to" part
          // but not anymore
          sync_update.cut_weight_to_after -= he_weight;
        }
      }

      // Get the new conductance objective
      HyperedgeWeight conductance_after = 
            ConductanceGlobalAttributedGains::gain(sync_update);
      tmp_scores[to] = conductance_after - conductance_now;
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
  // ! for the emulation of sync_update:
  // ! - to: the partition id of the node that is moved to
  // ! - cut_weight_from_after: the cut weight of the partition from after the move
  // ! - cut_weight_to_after: the cut weight of the partition to after the move
  // ! - volume_to_after: the used version of volume of the partition to after the move
  // ! Note: cut weights have to be adjusted later in the precomputeGains function
  template<typename PartitionedHypergraph>
  void setPreSyncUpdataData(const PartitionedHypergraph& phg,
                                const PartitionID to,
                                const bool original_stats,
                                SynchronizedEdgeUpdate& sync_update) {
    ASSERT(0 <= to && to <= _context.partition.k, "Invalid partition ID: " << V(to) << " k=" << _context.partition.k);
    ASSERT(sync_update.from != to, "Invalid partition ID: to=from=" << to << " k=" << _context.partition.k);
    ASSERT(0 <= sync_update.from && sync_update.from < _context.partition.k, "Invalid partition ID: from=" << sync_update.from << " k=" << _context.partition.k);
    sync_update.to = to;
    sync_update.cut_weight_from_after = phg.partCutWeight(sync_update.from); // will be adjusted
    sync_update.cut_weight_to_after = phg.partCutWeight(sync_update.to); // will be adjusted
    if (original_stats) {
      sync_update.volume_to_after = phg.partOriginalVolume(to) + sync_update.weighted_degree;
    } else {
      sync_update.volume_to_after = phg.partVolume(to) + sync_update.weighted_degree;
    }
  }

  using Base::_context;

  // ! Before gain computation, we construct a bitset that contains all
  // ! adjacent nodes of a block
  // [took from steiner_tree_gain_computation]
  tbb::enumerable_thread_specific<ds::Bitset> _local_adjacent_blocks;
  ds::Bitset _all_blocks;
};
 
}  // namespace mt_kahypar