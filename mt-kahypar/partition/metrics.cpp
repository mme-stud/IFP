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
// #include <tbb/parallel_reduce.h>

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/partition/mapping/target_graph.h"
#include "mt-kahypar/utils/exception.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"

namespace mt_kahypar::metrics {

namespace {

template<typename PartitionedHypergraph, Objective objective>
struct ObjectiveFunction { };

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::cut> {
  Gain operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return phg.connectivity(he) > 1 ? phg.edgeWeight(he) : 0;
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::km1> {
  Gain operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    return std::max(phg.connectivity(he) - 1, 0) * phg.edgeWeight(he);
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::soed> {
  Gain operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    const PartitionID connectivity = phg.connectivity(he);
    return connectivity > 1 ? connectivity * phg.edgeWeight(he) : 0;
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::steiner_tree> {
  Gain operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasTargetGraph());
    const TargetGraph* target_graph = phg.targetGraph();
    const Gain distance = target_graph->distance(phg.shallowCopyOfConnectivitySet(he));
    return distance * phg.edgeWeight(he);
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::conductance_local> {
  Gain operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasConductancePriorityQueue());
    ASSERT(phg.edgeIsEnabled(he), "Hyperedge " << he << " must be enabled to compute its contribution to the objective");

    const ds::ConductanceInfo top_conductance_info = phg.topPartConductanceInfo();
    const PartitionID top_part = top_conductance_info.partID;
    const HypernodeID pin_count_he_in_part = phg.pinCountInPart(he, top_part);
    ASSERT(0 <= top_part && top_part < phg.k());

    if (pin_count_he_in_part == 0 || pin_count_he_in_part == phg.edgeSize(he)) {
      // not a cutting edge => contributes by 0
      return 0;
    } else {
      // cutting edge => contributes by its weight / min(vol, total_vol - vol)
      const HypergraphVolume top_part_min_volume = top_conductance_info.fraction.getDenominator();
      // const HypergraphVolume total_volume_version = phg.conductancePriorityQueueUsesOriginalStats() ?
      //                                                  phg.originalTotalVolume() : 
      //                                                  phg.totalVolume();
      ASSERT(top_part_min_volume != 0);
      // ASSERT(total_volume_version != 0);


      // double_t scaled_contribution = 
      //    (static_cast<double_t>(phg.edgeWeight(he)) * static_cast<double_t>(total_volume_version)) /
      //    (static_cast<double_t>(top_part_min_volume) * static_cast<double_t>(phg.k()));
      double_t scaled_contribution = 
          (static_cast<double_t>(phg.edgeWeight(he)) * static_cast<double_t>(mt_kahypar::scaling_factor))
          / (static_cast<double_t>(top_part_min_volume));

      ASSERT(0 <= scaled_contribution);

      // HyperedgeWeight value_threshold = std::numeric_limits<HyperedgeWeight>::max();
      const Gain value_threshold = mt_kahypar::conductance_value_threshold;
      if (value_threshold < scaled_contribution) {
        LOG << "Scaled contribution of hyperedge " << he << " is too big: " << V(scaled_contribution) 
            << ". It is rounded to " << value_threshold;
        return value_threshold;
      }
      return static_cast<Gain>(scaled_contribution);
    }
    // TODO: several parts with the same conductance?
  }
};

template<typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::conductance_global> {
  Gain operator()(const PartitionedHypergraph& phg, const HyperedgeID& he) const {
    ASSERT(phg.hasConductancePriorityQueue());
    ASSERT(phg.edgeIsEnabled(he), "Hyperedge " << he << " must be enabled to compute its contribution to the objective");

    const ds::ConductanceInfo top_conductance_info = phg.topPartConductanceInfo();
    const PartitionID top_part = top_conductance_info.partID;
    const HypernodeID pin_count_he_in_part = phg.pinCountInPart(he, top_part);
    ASSERT(0 <= top_part && top_part < phg.k());

    if (pin_count_he_in_part == 0 || pin_count_he_in_part == phg.edgeSize(he)) {
      // not a cutting edge => contributes by 0
      return 0;
    } else {
      // cutting edge => contributes by its weight / min(vol, total_vol - vol)
      const HypergraphVolume top_part_min_volume = top_conductance_info.fraction.getDenominator();
      // const HypergraphVolume total_volume_version = phg.conductancePriorityQueueUsesOriginalStats() ?
      //                                                  phg.originalTotalVolume() : 
      //                                                  phg.totalVolume();
      ASSERT(top_part_min_volume != 0);
      // ASSERT(total_volume_version != 0);


      // double_t scaled_contribution = 
      //    (static_cast<double_t>(phg.edgeWeight(he)) * static_cast<double_t>(total_volume_version)) /
      //    (static_cast<double_t>(top_part_min_volume) * static_cast<double_t>(phg.k()));
      double_t scaled_contribution = 
          (static_cast<double_t>(phg.edgeWeight(he)) * static_cast<double_t>(mt_kahypar::scaling_factor))
          / (static_cast<double_t>(top_part_min_volume));

      ASSERT(0 <= scaled_contribution);

      // HyperedgeWeight value_threshold = std::numeric_limits<HyperedgeWeight>::max();
      const Gain value_threshold = mt_kahypar::conductance_value_threshold;
      if (value_threshold < scaled_contribution) {
        LOG << "Scaled contribution of hyperedge " << he << " is too big: " << V(scaled_contribution) 
            << ". It is rounded to " << value_threshold;
        return value_threshold;
      }
      return static_cast<Gain>(scaled_contribution);
    }
    // TODO: several parts with the same conductance?
  }
};

template <typename PartitionedHypergraph>
struct ObjectiveFunction<PartitionedHypergraph, Objective::aon_hypermodularity> {
  Gain operator()(const PartitionedHypergraph &phg,
                  const HyperedgeID &he) const {
    ASSERT(phg.hasAON());
    ASSERT(phg.hasOriginalEdgeSizes());
    std::size_t k = static_cast<std::size_t>(phg.originalEdgeSize(he));
    // LOG << "Hyperedge " << he << " with size " << k;
    if (k < 2) {
      // LOG << RED << "Size 1 edge " << he;
      return 0;                               // ignore size=1 edges
    }
    if (phg.connectivity(he) == 1) return 0;           // not cut

    long double beta_k = phg.beta(k);                        // β_k from finest level
    long double w_e    = static_cast<long double>(phg.edgeWeight(he));
    // LOG << " returns weights = " << w_e << " times " << beta_k;
    return static_cast<Gain>(w_e * beta_k );
  }
};

template <typename PartitionedHypergraph>
long double aonVolumeTerm(const PartitionedHypergraph& phg)
{
  const std::size_t dmax = static_cast<std::size_t>(phg.topLevelMaxEdgeSize());
  if (dmax < 2) return 0.0;

  long double term = 0.0;
  for (std::size_t d = 2; d <= dmax; ++d) {
    // LOG << "For edge size d = " << d << ": beta_d = " << phg.beta(d) << ", gamma_d = " << phg.gamma(d);
    long double inner = 0.0;
    for (auto c = 0; c < phg.k(); ++c) {
      inner += std::pow(phg.partOriginalVolume(c), static_cast<int>(d));
    }
    // term += phg.beta(d) * phg.gamma(d) * inner;
    term += phg.gamma(d) * inner;
  }
  return term;
}


template<typename PartitionedHypergraph>
Gain compute_conductance_objective(const PartitionedHypergraph& phg) {
  ASSERT( !PartitionedHypergraph::is_graph, "Conductance objective is not supported for graphs" );
  ASSERT(phg.hasConductancePriorityQueue());
  const ds::ConductanceInfo top_conductance_info = phg.topPartConductanceInfo();
  const HypergraphVolume top_part_cut_weight = top_conductance_info.fraction.getNumerator();
  const HypergraphVolume top_part_min_volume = top_conductance_info.fraction.getDenominator();

  if (top_part_min_volume < 1 || phg.k() < 2) {
    // only one block => bad partition
    LOG << "Top part min volume is too small: " << V(top_part_min_volume);
    return mt_kahypar::scaling_factor * 10;
  }
  // const HypergraphVolume total_volume_version = phg.conductancePriorityQueueUsesOriginalStats() ?
  //                                                      phg.originalTotalVolume() : 
  //                                                      phg.totalVolume();

  ASSERT(top_part_min_volume != 0);
  // ASSERT(total_volume_version != 0);
  ASSERT(top_part_cut_weight <= top_part_min_volume);


  // double_t scaled_conductance = 
  //    (static_cast<double_t>(top_part_cut_weight) * static_cast<double_t>(total_volume_version)) /
  //    (static_cast<double_t>(top_part_min_volume) * static_cast<double_t>(phg.k()));
  double_t scaled_conductance = static_cast<double_t>(top_part_cut_weight) 
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

template<Objective objective, typename PartitionedHypergraph>
Gain compute_objective_parallel(const PartitionedHypergraph& phg) {
  // Compute objective without for-loop through all hyperedges
  switch (objective) {
    case Objective::conductance_local:
    case Objective::conductance_global:
      return compute_conductance_objective(phg);
    default: break;
  }
  // Compute objective by iterating through all hyperedges
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  tbb::enumerable_thread_specific<Gain> obj(0);
  phg.doParallelForAllEdges([&](const HyperedgeID he) {
    obj.local() += func(phg, he);
  });
  
  if (objective != Objective::aon_hypermodularity) {
    return obj.combine(std::plus<>()) / (PartitionedHypergraph::is_graph ? 2 : 1);
  } else {
    Gain edge_sum = obj.combine(std::plus<>());
    Gain Q = edge_sum + aonVolumeTerm(phg);
    return Q;
  }
}

template<Objective objective, typename PartitionedHypergraph>
Gain compute_objective_sequentially(const PartitionedHypergraph& phg) {
  // Compute objective without for-loop through all hyperedges
  switch (objective) {
    case Objective::conductance_local:
    case Objective::conductance_global:
      return compute_conductance_objective(phg);
    default: break;
  }
  // Compute objective by iterating through all hyperedges
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  Gain obj = 0;
  for (const HyperedgeID& he : phg.edges()) {
    obj += func(phg, he);
  }

  
  if (objective != Objective::aon_hypermodularity) {
    return obj / (PartitionedHypergraph::is_graph ? 2 : 1);
  } else {
    Gain edge_sum = obj;
    Gain Q = edge_sum + aonVolumeTerm(phg);
    return Q;
  }
}

template<Objective objective, typename PartitionedHypergraph>
Gain contribution(const PartitionedHypergraph& phg, const HyperedgeID he) {
  ObjectiveFunction<PartitionedHypergraph, objective> func;
  return func(phg, he);
}

}

template<typename PartitionedHypergraph>
Gain quality(const PartitionedHypergraph& hg,
                        const Context& context,
                        const bool parallel) {
  return quality(hg, context.partition.objective, parallel);
}

template<typename PartitionedHypergraph>
Gain quality(const PartitionedHypergraph& hg,
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
    case Objective::aon_hypermodularity:
      return parallel ? compute_objective_parallel<Objective::aon_hypermodularity>(hg) :
        compute_objective_sequentially<Objective::aon_hypermodularity>(hg);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template<typename PartitionedHypergraph>
Gain contribution(const PartitionedHypergraph& hg,
                             const HyperedgeID he,
                             const Objective objective) {
  switch (objective) {
    case Objective::cut: return contribution<Objective::soed>(hg, he);
    case Objective::km1: return contribution<Objective::km1>(hg, he);
    case Objective::soed: return contribution<Objective::soed>(hg, he);
    case Objective::steiner_tree: return contribution<Objective::steiner_tree>(hg, he);
    case Objective::conductance_local: return contribution<Objective::conductance_local>(hg, he);
    case Objective::conductance_global: return contribution<Objective::conductance_global>(hg, he);
    case Objective::aon_hypermodularity: return contribution<Objective::aon_hypermodularity>(hg, he);
    default: throw InvalidParameterException("Unknown Objective");
  }
  return 0;
}

template<typename PartitionedHypergraph>
bool isBalanced(const PartitionedHypergraph& phg, const Context& context) {
  size_t num_empty_parts = 0;
  bool acceptable_part_weights = true;
  for (PartitionID i = 0; i < context.partition.k; ++i) {
    if (phg.partWeight(i) > context.partition.max_part_weights[i]) {
      acceptable_part_weights = false;
      break;
    }
    if (phg.partWeight(i) == 0) {
      num_empty_parts++;
    }
  }
  return (acceptable_part_weights && context.partition.preset_type == PresetType::large_k) ||
    (acceptable_part_weights && num_empty_parts <= phg.numRemovedHypernodes()) ||
    (context.partition.clustering && (phg.k() - num_empty_parts > 1));
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
    tbb::enumerable_thread_specific<Gain> approx_factor(0);
    hypergraph.doParallelForAllEdges([&](const HyperedgeID& he) {
      const size_t connectivity = hypergraph.connectivity(he);
      approx_factor.local() += connectivity <= context.mapping.max_steiner_tree_size ? 1 : 2;
    });
    return static_cast<double>(approx_factor.combine(std::plus<>())) / hypergraph.initialNumEdges();
  } else {
    return 1.0;
  }
}

template<typename PartitionedHypergraph>
double compute_double_conductance(const PartitionedHypergraph& phg) {
  ASSERT( !PartitionedHypergraph::is_graph, "Conductance objective is not supported for graphs" );
  ASSERT(phg.hasConductancePriorityQueue());
  const ds::ConductanceInfo top_conductance_info = phg.topPartConductanceInfo();
  const HypergraphVolume top_part_cut_weight = top_conductance_info.fraction.getNumerator();
  const HypergraphVolume top_part_min_volume = top_conductance_info.fraction.getDenominator();
  return static_cast<double_t>(top_part_cut_weight) /
         static_cast<double_t>(top_part_min_volume);
  /* 
    double conductance = tbb::parallel_reduce(
    tbb::blocked_range<PartitionID>(PartitionID(0), phg.k()),
      0, [&](const tbb::blocked_range<PartitionID>& range, double init) {
        if ( init >= 1 ) {
          return init;
        }
        double my_range_conductance = init;
        HypergraphVolume total_volume = phg.totalVolume();
        for (PartitionID p = range.begin(); p < range.end(); ++p) {
          HypergraphVolume part_volume = phg.partVolume(p);
          if (part_volume == 0) {
            continue;
          } else if (part_volume == total_volume) {
            my_range_conductance = 2;
            break;
          }
          HypergraphVolume part_min_volume = std::min(part_volume, total_volume - part_volume);
          HypergraphVolume part_cut_weight = phg.partCutWeight(p);    
          double part_conductance = static_cast<double>(part_cut_weight) /
                                    static_cast<double>(part_min_volume);
          my_range_conductance = std::max(my_range_conductance, part_conductance);
          ASSERT(1 >= my_range_conductance);
          if ( my_range_conductance == 1 ) {
            break;
          }
        }
        return my_range_conductance;
      }, [](const double lhs, const double rhs) {
            return std::max(lhs, rhs);
      });
  return conductance;
  */
}

template<typename PartitionedHypergraph>
long double compute_double_aon_hypermodularity(const PartitionedHypergraph& phg) {
  if (!phg.hasAON()) {
    return std::numeric_limits<long double>::quiet_NaN();
  }
  // Compute objective by iterating through all hyperedges
  auto double_edge_gain = [&](const PartitionedHypergraph& phg, const HyperedgeID he) 
    -> long double
    {
    std::size_t k = static_cast<std::size_t>(phg.originalEdgeSize(he));
    // LOG << "Hyperedge " << he << " with size " << k;
    if (k < 2) {
      return 0;                               // ignore size=1 edges
    }
    if (phg.connectivity(he) == 1) return 0;           // not cut

    long double beta_k = phg.beta(k);                        // β_k from finest level
    long double w_e    = static_cast<long double>(phg.edgeWeight(he));
    // LOG << " returns weights = " << w_e << " times " << beta_k;
    return w_e * beta_k;
  };

  long double edge_sum = 0.0;
  for (const HyperedgeID& he : phg.edges()) {
    edge_sum += double_edge_gain(phg, he);
  }
  long double Q = edge_sum + aonVolumeTerm(phg);
  return Q;

}

namespace {
#define OBJECTIVE_1(X) Gain quality(const X& hg, const Context& context, const bool parallel)
#define OBJECTIVE_2(X) Gain quality(const X& hg, const Objective objective, const bool parallel)
#define CONTRIBUTION(X) Gain contribution(const X& hg, const HyperedgeID he, const Objective objective)
#define IS_BALANCED(X) bool isBalanced(const X& phg, const Context& context)
#define IMBALANCE(X) double imbalance(const X& hypergraph, const Context& context)
#define APPROX_FACTOR(X) double approximationFactorForProcessMapping(const X& hypergraph, const Context& context)
#define CONDUCTANCE_DOUBLE(X) double compute_double_conductance(const X& phg)
#define AON_HYPERMODULARITY_DOUBLE(X) long double compute_double_aon_hypermodularity(const X& phg)
}

INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_1)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(OBJECTIVE_2)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONTRIBUTION)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IS_BALANCED)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(IMBALANCE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(APPROX_FACTOR)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(CONDUCTANCE_DOUBLE)
INSTANTIATE_FUNC_WITH_PARTITIONED_HG(AON_HYPERMODULARITY_DOUBLE)

} // namespace mt_kahypar::metrics