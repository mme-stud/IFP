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
 * The above copyright notice and this permission notice shall be included in
 *all copies or substantial portions of the Software.
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
#include "mt-kahypar/datastructures/sparse_map.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/refinement/gains/aon_hypermodularity/aon_hypermodularity_attributed_gains.h"
#include "mt-kahypar/partition/refinement/gains/gain_computation_base.h"
#include <unordered_set>
#include <vector>

namespace mt_kahypar {
class AONHypermodularityGainComputation
    : public GainComputationBase<AONHypermodularityGainComputation,
                                 AONHypermodularityAttributedGains> {
  using Base = GainComputationBase<AONHypermodularityGainComputation,
                                   AONHypermodularityAttributedGains>;
  using RatingMap = typename Base::RatingMap;

  static constexpr bool enable_heavy_assert = false;

public:
  AONHypermodularityGainComputation(const Context &context,
                                    bool disable_randomization = false)
      : Base(context, disable_randomization) {}

  // ! Precomputes the gain to all adjacent blocks.
  // ! Conceptually, we compute the gain of moving the node to an non-adjacent
  // block ! and the gain to all adjacent blocks assuming the node is in an
  // isolated block. ! The gain of that node to a block to can then be computed
  // by ! 'isolated_block_gain - tmp_scores[to]' (see gain(...))
  template <typename PartitionedHypergraph>
  void precomputeGains(const PartitionedHypergraph &phg, const HypernodeID hn,
                       RatingMap &tmp_scores, Gain &isolated_block_gain,
                       const bool) {
    ASSERT(tmp_scores.size() == 0, "Rating map not empty");
    const PartitionID from = phg.partID(hn);
    const double dv = static_cast<double>(phg.nodeOriginalWeightedDegree(hn));

    /* ====================================================================== *
     *  1. gather all adjacent blocks of hn
     * ====================================================================== */

    std::unordered_set<PartitionID> neighbs;
    for (HyperedgeID he : phg.incidentEdges(hn)) {
      for (auto clus : phg.connectivitySet(he)) {
        neighbs.insert(clus);
      }
    }

    /* ====================================================================== *
     *  2. compute dCut for all candidate blocks
     * ====================================================================== */

    vec<double> cut_scores(_context.partition.k, 0.0); // 1 slot / block

    for (const HyperedgeID he : phg.incidentEdges(hn)) {
      const std::size_t d = static_cast<std::size_t>(phg.originalEdgeSize(he));
      // HypernodeID computed_size = 0;
      // for (HypernodeID pin : phg.pins(he)) {
      //   computed_size += phg.nodeWeight(pin);
      // }
      // if(computed_size != phg.edgeStrength(he)) {
      //   LOG << RED << "Inv";
      // }
      // LOG << "Hyperedge " << he << " has size " << d << " and strength " << phg.edgeStrength(he);
      
      
      if (d < 2)
        continue;

      const double w_beta = phg.edgeWeight(he) * phg.beta(d);
      const PartitionID conn = phg.connectivity(he);
      const HypernodeID pins_from = phg.pinCountInPart(he, from);

      if (conn == 1 && d > 1) {
        // edge is internal → becomes cut no matter *where* we move
        for (const PartitionID blk : neighbs)
          if (blk != from)
            cut_scores[blk] += w_beta;
      } else if (conn == 2 && pins_from == 1) {
        // edge touches exactly two blocks and *hn* is the only pin in <from>
        // -> moving hn to the *other* block makes the edge internal
        for (const PartitionID blk : phg.connectivitySet(he))
          if (blk != from)
            cut_scores[blk] -= w_beta; // blk == "other" (= candidate)
      }
    }

    /* ------------------------------------------------------------
     * 3. iterate once over every candidate block
     * ---------------------------------------------------------- */

    /* volumes (needed for ΔVol) */
    const std::size_t d_max = std::max({
        static_cast<std::size_t>(phg.topLevelMaxEdgeSize()),
        phg.betaVector().size() - 1,
        phg.gammaVector().size() - 1
    });
    const double vol_from = static_cast<double>(phg.partOriginalVolume(from));

    for (PartitionID to : neighbs) {
      if (from == to)
        continue;

      /* ---- ΔVol -------------------------------------------------- */
      double dVol = 0.0;
      for (std::size_t d = 2; d <= d_max; ++d) {
        double g = phg.gamma(d);
        // double b = phg.beta(d);
        if (g == 0.0)
          continue;
        double vol_to = static_cast<double>(phg.partOriginalVolume(to));
        double delta = + std::pow(vol_from - dv, (int)d)
                       + std::pow(vol_to + dv, (int)d)
                       - std::pow(vol_from, (int)d)
                       - std::pow(vol_to, (int)d);
        // dVol += (b * g) * delta; // <- sign fixed here
        dVol += g * delta; // <- sign fixed here
      }

      /* ---- ΔCut -------------------------------------------------- */
      // double dCut = 0.0;

      // for (HyperedgeID he : phg.incidentEdges(hn)) {
      //   std::size_t d = static_cast<std::size_t>(phg.edgeSize(he));
      //   if (d < 2)
      //     continue;

      //   const double wbeta = phg.edgeWeight(he) * phg.beta(d);

      //   /* check if e\{hn} is homogeneous and record its cluster p1 */
      //   PartitionID p1 = kInvalidPartition;
      //   bool multi = false;
      //   for (HypernodeID pin : phg.pins(he)) {
      //     if (pin == hn)
      //       continue;
      //     PartitionID pid = phg.partID(pin);
      //     if (p1 == kInvalidPartition)
      //       p1 = pid;
      //     else if (pid != p1) {
      //       multi = true;
      //       break;
      //     }
      //   }
      //   if (multi)
      //     continue; // already cut → Δ = 0

      //   if (p1 == from) { // edge internal now
      //     if (to != p1)
      //       dCut += wbeta; // moving cuts it
      //   } else {           // edge already cut two-way
      //     if (to == p1)
      //       dCut -= wbeta; // moving un-cuts it
      //   }
      // }

      tmp_scores[to] = -(cut_scores[to] + dVol); // store total change
    }

    isolated_block_gain = 0; // baseline fixed to zero
  }

  Gain gain(const Gain to_score, const Gain isolated_block_gain) {
    return isolated_block_gain - to_score;
  }

  void changeNumberOfBlocksImpl(const PartitionID) {
    // Do nothing
  }
};
} // namespace mt_kahypar