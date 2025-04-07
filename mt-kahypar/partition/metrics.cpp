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

#include "mt-kahypar/partition/metrics.h"

#include <cmath>
#include <algorithm>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/exception.h"
// include "mt-kahypar/datastructures/conductance_pq.h"

namespace mt_kahypar::metrics {

namespace {

template<typename PartitionedHypergraph, Objective objective>
struct ObjectiveFunction { };

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::cut> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return phg.connectivity(he) > 1 ? phg.edgeWeight(he) : 0;
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::km1> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return std::max(phg.connectivity(he) - 1, 0) * phg.edgeWeight(he);
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::soed> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    const PartitionID connectivity = phg.connectivity(he);
    return connectivity > 1 ? connectivity * phg.edgeWeight(he) : 0;
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::steiner_tree> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasTargetGraph());
    const TargetGraph* target_graph = phg.targetGraph();
    const HyperedgeWeight distance = target_graph->distance(phg.shallowCopyOfConnectivitySet(he));
    return distance * phg.edgeWeight(he);
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::conductance_local> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasConductancePriorityQueue());
    ASSERT(phg.edgeIsEnabled(he), "Hyperedge " << he << " must be enabled to compute its contribution to the objective");

    const ConductanceInfo top_conductance_info = phg.topPartConductanceInfo();
    const PartitionID top_part = top_conductance_info.partID;
    const HypernodeID pin_count_he_in_part = phg.pinCountInPart(he, top_part);
    ASSERT(0 <= top_part && top_part < phg.k());

    if (pin_count_he_in_part == 0 || pin_count_he_in_part == phg.edgeSize(he)) {
      // not a cutting edge => contributes by 0
      return 0;
    } else {
      // cutting edge => contributes by its weight / min(vol, total_vol - vol)
      const HypergraphVolume top_part_min_volume = top_conductance_info.fraction.getDenominator();
      const HypergraphVolume total_volume_version = phg.conductancePriorityQueueUsesOriginalStats() ?
                                                        phg.originalTotalVolume() : 
                                                        phg.totalVolume();
      ASSERT(top_part_min_volume != 0);
      ASSERT(total_volume_version != 0);

      double_t conductance_contribution = static_cast<double_t>(phg.edgeWeight(he)) /
                                          static_cast<double_t>(top_part_min_volume);
      double_t current_multiplier = static_cast<double_t>(phg.totalVolume()) / 
                                    static_cast<double_t>(phg.k());
      // double_t scaled_contribution = conductance_contribution * current_multiplier;
      double_t scaled_contribution = 
          (static_cast<double_t>(phg.edgeWeight(he)) * static_cast<double_t>(total_volume_version)) /
          (static_cast<double_t>(top_part_min_volume) * static_cast<double_t>(phg.k()));
      
      ASSERT(0 <= scaled_contribution);

      HyperedgeWeight value_threshold = std::numeric_limits<HyperedgeWeight>::max();
      if (value_threshold < scaled_contribution) {
        LOG << "Scaled of hyperedge " << he << " is too big: " << V(scaled_contribution) 
            << ". It is rounded to " << value_threshold;
        return value_threshold;
      }
      return static_cast<HyperedgeWeight>(scaled_contribution);
    }
    // TODO: several parts with the same conductance?
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::conductance_global> {
  HyperedgeWeight operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasConductancePriorityQueue());
    ASSERT(phg.edgeIsEnabled(he), "Hyperedge " << he << " must be enabled to compute its contribution to the objective");

    const ConductanceInfo top_conductance_info = phg.topPartConductanceInfo();
    const PartitionID top_part = top_conductance_info.partID;
    const HypernodeID pin_count_he_in_part = phg.pinCountInPart(he, top_part);
    ASSERT(0 <= top_part && top_part < phg.k());

    if (pin_count_he_in_part == 0 || pin_count_he_in_part == phg.edgeSize(he)) {
      // not a cutting edge => contributes by 0
      return 0;
    } else {
      // cutting edge => contributes by its weight / min(vol, total_vol - vol)
      const HypergraphVolume top_part_min_volume = top_conductance_info.fraction.getDenominator();
      const HypergraphVolume total_volume_version = phg.conductancePriorityQueueUsesOriginalStats() ?
                                                        phg.originalTotalVolume() : 
                                                        phg.totalVolume();
      ASSERT(top_part_min_volume != 0);
      ASSERT(total_volume_version != 0);

      double_t conductance_contribution = static_cast<double_t>(phg.edgeWeight(he)) /
                                          static_cast<double_t>(top_part_min_volume);
      double_t current_multiplier = static_cast<double_t>(phg.totalVolume()) / 
                                    static_cast<double_t>(phg.k());
      // double_t scaled_contribution = conductance_contribution * current_multiplier;
      double_t scaled_contribution = 
          (static_cast<double_t>(phg.edgeWeight(he)) * static_cast<double_t>(total_volume_version)) /
          (static_cast<double_t>(top_part_min_volume) * static_cast<double_t>(phg.k()));
      
      ASSERT(0 <= scaled_contribution);

      HyperedgeWeight value_threshold = std::numeric_limits<HyperedgeWeight>::max();
      if (value_threshold < scaled_contribution) {
        LOG << "Scaled of hyperedge " << he << " is too big: " << V(scaled_contribution) 
            << ". It is rounded to " << value_threshold;
        return value_threshold;
      }
      return static_cast<HyperedgeWeight>(scaled_contribution);
    }
    // TODO: several parts with the same conductance?
  }
};

template<Objective objective, typename PartitionedHypergraph>
HyperedgeWeight compute_objective_parallel(const PartitionedHypergraph& phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  tbb::enumerable_thread_specific<HyperedgeWeight> obj(0);
  phg.doParallelForAllEdges([&](const HyperedgeID he) {
    obj.local() += func(phg, he);
  });
  return obj.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template<Objective objective, typename PartitionedHypergraph>
HyperedgeWeight compute_objective_sequentially(const PartitionedHypergraph& phg) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  HyperedgeWeight obj = 0;
  for (const HyperedgeID& he : phg.edges()) {
    obj += func(phg, he);
  }
  return obj / (PartitionedHypergraph::is_graph ? 2 : 1);
}

template<Objective objective, typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& phg, const HyperedgeID he) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  return func(phg, he);
}

}

template<typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Context& context,
                        const bool parallel) {
  return quality(hg, context.partition.objective, parallel);
}

template<typename PartitionedHypergraph>
HyperedgeWeight quality(const PartitionedHypergraph& hg,
                        const Objective objective,
                        const bool parallel) {
  switch (objective) {
    case Objective::cut:
      return parallel ? compute_objective_parallel<Objective::cut>(hg) :
        compute_objective_sequentially<Objective::cut>(hg);
    case Objective::km1:
      return parallel ? compute_objective_parallel<Objective::km1>(hg) :
        compute_objective_sequentially<Objective::km1>(hg);
    case Objective::soed:
      return parallel ? compute_objective_parallel<Objective::soed>(hg) :
        compute_objective_sequentially<Objective::soed>(hg);
    case Objective::steiner_tree:
      return parallel ? compute_objective_parallel<Objective::steiner_tree>(hg) :
        compute_objective_sequentially<Objective::steiner_tree>(hg);
    case Objective::conductance_local:
      return parallel ? compute_objective_parallel<Objective::conductance_local>(hg) :
        compute_objective_sequentially<Objective::conductance_local>(hg);
    case Objective::conductance_global:
      return parallel ? compute_objective_parallel<Objective::conductance_global>(hg) :
        compute_objective_sequentially<Objective::conductance_global>(hg);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template<typename PartitionedHypergraph>
HyperedgeWeight contribution(const PartitionedHypergraph& hg,
                             const HyperedgeID he,
                             const Objective objective) {
  switch (objective) {
    case Objective::cut: return contribution<Objective::soed>(hg, he);
    case Objective::km1: return contribution<Objective::km1>(hg, he);
    case Objective::soed: return contribution<Objective::soed>(hg, he);
    case Objective::steiner_tree: return contribution<Objective::steiner_tree>(hg, he);
    case Objective::conductance_local: return contribution<Objective::conductance_local>(hg, he);
    case Objective::conductance_global: return contribution<Objective::conductance_global>(hg, he);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template<typename PartitionedHypergraph>
bool isBalanced(const PartitionedHypergraph& phg, const Context& context) {
  size_t num_empty_parts = 0;
  for (PartitionID i = 0; i < context.partition.k; ++i) {
    if (phg.partWeight(i) > context.partition.max_part_weights[i]) {
      return false;
    }
    if (phg.partWeight(i) == 0) {
      num_empty_parts++;
    }
  }
  return context.partition.preset_type == PresetType::large_k ||
    num_empty_parts <= phg.numRemovedHypernodes();
}

template<typename PartitionedHypergraph>
double imbalance(const PartitionedHypergraph& hypergraph, const Context& context) {
  ASSERT(context.partition.perfect_balance_part_weights.size() == (size_t)context.partition.k);

  double max_balance = (hypergraph.partWeight(0) /
                        static_cast<double>(context.partition.perfect_balance_part_weights[0]));

  for (PartitionID i = 1; i < context.partition.k; ++i) {
    const double balance_i =
            (hypergraph.partWeight(i) /
              static_cast<double>(context.partition.perfect_balance_part_weights[i]));
    max_balance = std::max(max_balance, balance_i);
  }

  return max_balance - 1.0;
}

template<typename PartitionedHypergraph>
double approximationFactorForProcessMapping(const PartitionedHypergraph& hypergraph, const Context& context) {
  if ( !PartitionedHypergraph::is_graph ) {
    tbb::enumerable_thread_specific<HyperedgeWeight> approx_factor(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      const size_t connectivity = hypergraph.connectivity(he);
      approx_factor.local() += connectivity <= context.mapping.max_steiner_tree_size ? 1 : 2;
    });
    return static_cast<double>(approx_factor.combine(std::plus<>())) / hypergraph.initialNumEdges();
  } else {
    return 1.0;
  }
}

namespace {
#define OBJECTIVE_1(X) HyperedgeWeight quality(const X& hg, const Context& context, const bool parallel)
#define OBJECTIVE_2(X) HyperedgeWeight quality(const X& hg, const Objective objective, const bool parallel)
#define CONTRIBUTION(X) HyperedgeWeight contribution(const X& hg, const HyperedgeID he, const Objective objective)
#define IS_BALANCED(X) bool isBalanced(const X& phg, const Context& context)
#define IMBALANCE(X) double imbalance(const X& hypergraph, const Context& context)
#define APPROX_FACTOR(X) double approximationFactorForProcessMapping(const X& hypergraph, const Context& context)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_1)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_2)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONTRIBUTION)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IS_BALANCED)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IMBALANCE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(APPROX_FACTOR)

} // namespace mt_kahypar::metrics