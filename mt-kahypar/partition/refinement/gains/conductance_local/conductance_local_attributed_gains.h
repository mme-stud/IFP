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

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/conductance_pq.h"
#include "mt-kahypar/datastructures/nonnegative_fraction.h"

namespace mt_kahypar {

/**
 * After moving a node, we perform a synchronized update of the pin count values
 * for each incident hyperedge of the node based on which we then compute an
 * attributed gain value.
 */
struct ConductanceLocalAttributedGains {
  // sync_update must contain the following values:
  // - cut_weight_from_before: the cut weight of the partition from before the move
  // - cut_weight_from_after: the cut weight of the partition from after the move
  // - cut_weight_to_before: the cut weight of the partition to before the move
  // - cut_weight_to_after: the cut weight of the partition to after the move
  // - volume_from_after: the used version of volume of the partition from after the move
  // - volume_to_after: the used version of volume of the partition to after the move
  // - total_volume: the used version of total volume of the hypergraph
  // - weighted_degree: the used version of weighted degree of the node that is moved
  // - k: the number of partitions
  static Gain gain(const SynchronizedEdgeUpdate& sync_update) {
    // ASSERT(SyncUpdatePreferences::collective_sync_updates_in_phg, 
    //  "Synchronized gain updates should be enabled for Conductance Attribted Gains");
    // Note: can't check if collective sync_updates are enabled, 
    //       but PartitionedHypergraph::topThreePartConductanceInfos() asserts it

    // Find old local max
    ds::ConductanceFraction old_fraction_from(
      sync_update.cut_weight_from_before,
      std::min(sync_update.volume_from_after + sync_update.weighted_degree, 
               sync_update.total_volume - (sync_update.volume_from_after + sync_update.weighted_degree))
    );
    ds::ConductanceFraction old_fraction_to(
      sync_update.cut_weight_to_before,
      std::min(sync_update.volume_to_after - sync_update.weighted_degree, 
               sync_update.total_volume - (sync_update.volume_to_after - sync_update.weighted_degree))
    );
    ds::ConductanceFraction old_local_top_fraction = old_fraction_from;
    if (old_fraction_to > old_local_top_fraction) {
      old_local_top_fraction = old_fraction_to;
    }

    // Find new local max
    ds::ConductanceFraction new_fraction_from(
      sync_update.cut_weight_from_after,
      std::min(sync_update.volume_from_after, sync_update.total_volume - sync_update.volume_from_after)
    );
    ds::ConductanceFraction new_fraction_to(
      sync_update.cut_weight_to_after,
      std::min(sync_update.volume_to_after, sync_update.total_volume - sync_update.volume_to_after)
    );
    // The old top conductance fraction
    ds::ConductanceFraction new_local_top_fraction = new_fraction_from;
    if (new_fraction_to > new_local_top_fraction) {
      new_local_top_fraction = new_fraction_to;
    }

    Gain new_local_conductance = compute_conductance_objective(sync_update.total_volume,
                                                                          new_local_top_fraction,
                                                                          sync_update.k);
    Gain old_local_conductance = compute_conductance_objective(sync_update.total_volume,
                                                                          old_local_top_fraction,
                                                                          sync_update.k);
    return new_local_conductance - old_local_conductance; // gain is positive if conductance increases
  }

  static Gain compute_conductance_objective(const HypergraphVolume& total_volume_version,
                                                       const ds::ConductanceFraction& fraction, 
                                                       const PartitionID& k) {
    unused(k);
    unused(total_volume_version);
    const HypergraphVolume& top_part_cut_weight = fraction.getNumerator();
    const HypergraphVolume& top_part_min_volume = fraction.getDenominator();
    if ( top_part_min_volume == 0 ) {
      // return std::numeric_limits<HyperedgeWeight>::max();
      return mt_kahypar::scaling_factor * 10;
    }
    ASSERT(top_part_min_volume != 0);
  
    //if ( top_part_cut_weight > top_part_min_volume ) {
    //  LOG << "Weights are wrong: " << V(top_part_cut_weight) << " > " << V(top_part_min_volume);
    //} else {
    //  LOG << "Weights are right";
    //}
    if ( top_part_cut_weight <= top_part_min_volume ) {
      Gain scaled_conductance = static_cast<double_t>(top_part_cut_weight) 
                                  / static_cast<double_t>(top_part_min_volume)
                                  * static_cast<double_t>(mt_kahypar::scaling_factor);
      ASSERT(0 <= scaled_conductance);
      // HyperedgeWeight value_threshold = std::numeric_limits<HyperedgeWeight>::max();
      const Gain value_threshold = mt_kahypar::conductance_value_threshold;
      if (value_threshold < scaled_conductance) {
        LOG << "Scaled conductance is too big: " << V(scaled_conductance) 
            << ". It is rounded to " << value_threshold;
        return value_threshold;
      }
      return static_cast<Gain>(scaled_conductance);
    }
    // weights are wrong, but conductnace could  therefore be high
    return mt_kahypar::scaling_factor;
  }

};

}  // namespace mt_kahypar
 