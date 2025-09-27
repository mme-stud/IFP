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

// Adil's code

#pragma once

#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar {
/**
 * After moving a node, we perform a synchronized update of the pin count values
 * for each incident hyperedge of the node based on which we then compute an
 * attributed gain value.
 */
struct AONHypermodularityAttributedGains {
  static Gain gain(const SynchronizedEdgeUpdate &sync_update) {
    long double dVol = volumeDelta(sync_update) / sync_update.hn_degree;
    // LOG << "Moving node " << sync_update.hn << " from " << sync_update.from
    // << " to " << sync_update.to << "; dvol = " << volumeDelta(sync_update);
    if (sync_update.edge_size < 2) // not a cutting edge
      return dVol;
    const vec<AONCoefficient> &beta = *sync_update.beta_vec;
    const long double wbeta = (
        sync_update.original_edge_size < beta.size() ?
          sync_update.edge_weight * beta[sync_update.original_edge_size] :
          0.0);

    //  (a) edge was internal → becomes cut  (+)
    if ((sync_update.pin_count_in_from_part_after + 1) ==
        sync_update.edge_size) { // cut after
      return wbeta + dVol;
    }

    //  (b) edge was 2-way cut → becomes internal again (–)
    if (sync_update.pin_count_in_to_part_after ==
        sync_update.edge_size) { // internal after
      return -wbeta + dVol;
    }

    return dVol;
  }

  static long double volumeDelta(const SynchronizedEdgeUpdate &sync_update) {
    long double dVol = 0;
    // const vec<AONCoefficient> &beta = *sync_update.beta_vec;
    const vec<AONCoefficient> &gamma = *sync_update.gamma_vec;
    const long double v_from = static_cast<long double>(
        sync_update.original_volume_from_after + sync_update.original_weighted_degree); // volume *before*
    const long double v_to = static_cast<long double>(
        sync_update.original_volume_to_after - sync_update.original_weighted_degree); // volume *before*
    const long double dv =
        static_cast<long double>(sync_update.original_weighted_degree); // (= node volume)

    std::size_t max_d = std::min({ // beta.size() - 1,
                                   gamma.size() - 1,
                                   static_cast<std::size_t>(sync_update.max_edge_size)});
    for (std::size_t d = 2; d <= max_d; ++d) {

      // const long double b = beta[d];
      const long double g = gamma[d];
      if (g == 0.0)
        continue;

      const long double delta = + std::pow(v_from - dv, int(d))
                           + std::pow(v_to + dv, int(d))
                           - std::pow(v_from, int(d))
                           - std::pow(v_to, int(d));
      // dVol += b * g * delta;
      dVol += g * delta;
    }
    return dVol;
  }
};
} // namespace mt_kahypar