/*******************************************************************************
 * MIT License
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

#include "mt-kahypar/partition/initial_partitioning/aon_hypermodularity_initial_partitioner.h"
#include "mt-kahypar/datastructures/sequential_ip_hypergraph.h"
#include "mt-kahypar/datastructures/static_hypergraph.h"

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"

#include <tbb/task_arena.h>
 
/// [debug] #include "mt-kahypar/IFP/mt-kahypar/io/hypergraph_io.h"
 
namespace mt_kahypar {

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::partitionImpl(
    const HypernodeID edgeSizeThreshold, 
    const long long maxNumIter, 
    const AONCoefficient eps, 
    const AONCoefficient clusterPenalty,
    const bool randomize,
    bool useOriginalEdgeSizes
    ) {
  /// [debug] std::cout << "partitionImpl()" << std::endl;
  // if num. nodes = k, assign each node to its own block
  // otherwise produce same result as random IP (maybe change later)
  if ( _ip_data.should_initial_partitioner_run(_ipName) ) {
  tbb::task_arena limited_arena(1);
  limited_arena.execute([&]{
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
    if (_ip_data._is_local_hg_initialized.local() != -1) {
      LOG << RED << "AON IP [" << V(_ipName) << V(_tag) << "]: Local hypergraph already initialized! " << _ip_data._is_local_hg_initialized.local() << END;
    }
    _ip_data._is_local_hg_initialized.local() = _tag;
    LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: Starting";

    if (hg.hasFixedVertices() || !hg.hasAON() || !hg.is_static_hypergraph || hg.is_graph) {
      // fixed vertices are currently not supported in AON-Hypermodularity IP
      randomPartition(hg);
    } else {
      AON_HMLL(hg, 
               edgeSizeThreshold,
               maxNumIter,
               eps,
               clusterPenalty,
               randomize,
               useOriginalEdgeSizes);
    }
    // =============== General final steps of IP =============

    hg.needsConductancePriorityQueue();

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(_ipName, _rng, _tag, time);
    if (_ip_data._is_local_hg_initialized.local() != _tag) 
      LOG << RED << "AON IP [" << V(_ipName) << V(_tag) << "]: Local hypergraph not initialized after commit! " << _ip_data._is_local_hg_initialized.local() << END;
    _ip_data._is_local_hg_initialized.local() = -1;
    LOG << "AON IP [" << V(_ipName) << V(_tag) << "] finished in " << time << " seconds.";
  }
  );
  }
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::AON_HMLL(PartitionedHypergraph& hg,
                                                                const HypernodeID edgeSizeThreshold, 
                                                                const long long maxNumIter, 
                                                                const AONCoefficient eps, 
                                                                const AONCoefficient clusterPenalty,
                                                                const bool randomize,
                                                                bool useOriginalEdgeSizes) {
  // coarsest hypergraph - sequential version
  SequentialIPHypergraph H_coarsest(hg.hypergraph());

  if (! useOriginalEdgeSizes || ! H_coarsest.hasOriginalEdgeSizes()) {
    useOriginalEdgeSizes = false;
    if ( ! H_coarsest.hasOriginalEdgeSizes() && useOriginalEdgeSizes)
      if (_context.partition.verbose_output)
        LOG << "No snapshot of original edge size found: using edge sizes of the current hypergraph.";
    // save current edge sizes, weighted degrees and total volume
    H_coarsest.snapshotOriginalEdgeSizes();
    H_coarsest.snapshotOriginalWeightedDegreesAndTotalVolume();
  }
  if (_context.partition.verbose_output) {
    if (useOriginalEdgeSizes)
      LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: Using original edge sizes: # he " << H_coarsest.initialNumEdges() << "# hn " << H_coarsest.initialNumNodes();
    else
      LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: Using current edge sizes. # he " << H_coarsest.initialNumEdges() << "# hn " << H_coarsest.initialNumNodes();
  }
  H_coarsest.enableSinglePinNetsRemoval(); // single pin nets are never cutting
  H_coarsest.useOriginalSizeInParallelNetsDetection(true); // otherwise gain is incorrect
  
  // =====================================================
  //          1. Singleton initial partitioning
  // =====================================================
  
    // current communities of hg: z: node -> community
    std::vector<PartitionID> z(hg.initialNumNodes(), kInvalidPartition);
  for (const HypernodeID &hn : H_coarsest.nodes()) {
    z[hn] = hn;
    // H_coarsest.setCommunityID(hn, hn); // set in constructor
  }

  // =====================================================
  //          2. AllOrNothingHMLL: Louvain Cycle
  // =====================================================
  bool z_changed = false;
  AONCoefficient total_gain = 0.0;
  long long counter = 0;
  do {
    counter++;
    /// [debug] std::cout << "Outer Iteration: " << counter << std::endl;

    /** ------------------ Louvain Step: ------------------
     *  - Nodes are moved to neighboring partitions as
     *    long as it improves the modularity gain;
     */
    AONCoefficient new_gain = louvainStep(H_coarsest,
                                          edgeSizeThreshold, 
                                          maxNumIter, eps, clusterPenalty, randomize);
    if (_context.partition.verbose_output)
      LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: step " << counter << " gain = " << new_gain;
    total_gain += new_gain;
    z_changed = (new_gain > eps);
    /// [debug] std::cout << "Outer Iteration: louvainStep(..) finished " << counter << std::endl;

    /** --------------------- Expand: ---------------------
     *  - z (communities for H) is updated via contract();
     */
    /** -------------------- Collapse: --------------------
     *  - The partitions of H_coarsest are collapsed into
     *  nodes;
     */
    expandAndCollapse(H_coarsest, z);
    /// [debug] std::cout << "Outer Iteration: expandAndCollapse(..) finished " << counter << std::endl;
  } while (z_changed);

  if (_context.partition.verbose_output)
    LOG << "AON IP [" << V(_ipName) << V(_tag) << "] finished Louvain with total gain " << total_gain 
        << " and " << H_coarsest.k() << " clusters.";

  // =====================================================
  //             3. Finalize Partitioning
  // =====================================================

  // HypernodeID new_k = H_coarsest.k();
  // hg.setK(new_k); 
  for (const HypernodeID &hn : hg.nodes()) {
    PartitionID partition = z[hn];
    ASSERT(partition != kInvalidPartition,
        "AONHypermodularityInitialPartitioner::partitionImpl: "
        "node " << hn << " is not assigned to a new partition");
    ASSERT(partition < H_coarsest.k(), "AONHypermodularityInitialPartitioner::partitionImpl: "
        "node " << hn << " is assigned to an invalid partition: k = " << H_coarsest.k());
    hg.setOnlyNodePart(hn, partition);
  }
  hg.initializePartition(false /* parallel */);
}

template<typename TypeTraits>
AONCoefficient AONHypermodularityInitialPartitioner<TypeTraits>::louvainStep(
            SequentialIPHypergraph& H_coarsest,
            const HypernodeID edgeSizeThreshold, const long long maxNumIter, const AONCoefficient eps, const AONCoefficient clusterPenalty, const bool randomize) {
  // precompute neighboring nodes
  HypernodeID numNodes = H_coarsest.initialNumNodes();
  ASSERT(H_coarsest.k() == static_cast<PartitionID>(numNodes), 
        "AONHypermodularityInitialPartitioner::louvainStep: "
        "no singleton partition at the start: "
        "k = " << H_coarsest.k() << ", numNodes = " << numNodes);
  std::vector<std::vector<HypernodeID>> neighbors(numNodes, std::vector<HypernodeID>(0));
  std::vector<bool> visited(numNodes, false);
  for (const HypernodeID &i : H_coarsest.nodes()) {
    visited.assign(numNodes, false);
    for (const HyperedgeID &he : H_coarsest.incidentEdges(i)) {
      for (const HypernodeID &hn : H_coarsest.pins(he)) {
        ASSERT(H_coarsest.partID(hn) == static_cast<PartitionID>(hn), 
              "AONHypermodularityInitialPartitioner::louvainStep: "
              "node " << hn << " is not assigned to its singleton partition: "
              "partID(hn) = " << H_coarsest.partID(hn));
        visited[hn] = true;
      }
    }
    visited[i] = false;
    for (const HypernodeID& hn : H_coarsest.nodes()) {
      if (visited[hn])
        neighbors[i].push_back(hn);
    }

    ASSERT(static_cast<PartitionID>(neighbors[i].size()) < H_coarsest.k(), 
           "AONHypermodularityInitialPartitioner::louvainStep: "
           "node " << i << " has more neighboring nodes than partitions");
  }

  AONCoefficient total_gain = 0.0;
  bool improving = true;
  long long iter = 0;
  PartitionID k_now = H_coarsest.k(); // tracks current number of clusters
  if (randomize) {
    std::vector<HypernodeID> nodes(numNodes, 0);
    for (HypernodeID i = 0; i < numNodes; ++i) {
      nodes[i] = i;
    }
    std::shuffle(nodes.begin(), nodes.end(), _rng);
    while (improving && (iter++ < maxNumIter)) {
      AONCoefficient gain = 0.0;
      improving = false;
      for (const HypernodeID &i : nodes) {
        if ( ! H_coarsest.nodeIsEnabled(i) ) continue;
        gain += louvainStepForANode(i, neighbors[i], visited, 
                                       k_now /* reference! */,
                                       H_coarsest,
                                       edgeSizeThreshold, 
                                       eps, clusterPenalty, randomize);
      }
      total_gain += gain;
      improving = (gain > eps);
    }
  } else {
    while (improving && (iter++ < maxNumIter)) {
      AONCoefficient gain = 0.0;
      improving = false;
      /// [debug]   std::cout << "Louvain: round " << iter << std::endl;
      for (const HypernodeID &i : H_coarsest.nodes()) {
        gain += louvainStepForANode(i, neighbors[i], visited, 
                                    k_now /* reference! */, 
                                    H_coarsest,
                                    edgeSizeThreshold, 
                                    eps, clusterPenalty, randomize);
      }
      total_gain += gain;
      improving = (gain > eps);
    }
  }
  return total_gain;
}

template<typename TypeTraits>
AONCoefficient AONHypermodularityInitialPartitioner<TypeTraits>::louvainStepForANode(
        const HypernodeID& i, const std::vector<HypernodeID>& neighbors_i, std::vector<bool>& visitedParts, 
        PartitionID& k_now,
        SequentialIPHypergraph& H_coarsest,
        const HypernodeID edgeSizeThreshold, const AONCoefficient eps, const AONCoefficient clusterPenalty, const bool randomize) {
  unused(randomize);
  /// [debug] if (i % 1000 == 0)
  /// [debug] std::cout << "Louvain: node " << i << std::endl;
  PartitionID part_i = H_coarsest.partID(i);

  /// Check all neighboring partitions to find the best gain
  visitedParts.assign(H_coarsest.k(), false);
  visitedParts[part_i] = true; // mark current partition as visited

  AONCoefficient cluster_gain = 0.0;
  if (clusterPenalty > 0.0 && H_coarsest.partSize(part_i) == 1) {
    cluster_gain = clusterPenalty * (std::log(k_now) - std::log(k_now - 1));
  }

  AONCoefficient best_gain = 0.0;
  PartitionID best_partition = part_i;
  // for (const HyperedgeID &he : H_coarsest.incidentEdges(i)) {
  //   for (const PartitionID &A : H_coarsest.connectivitySet(he)) {
  // bool improving = false;
  for (const HypernodeID &neighbor : neighbors_i) {
    PartitionID A = H_coarsest.partID(neighbor);
    if (visitedParts[A]) continue;
    visitedParts[A] = true;
    AONCoefficient gain = QAONGain(H_coarsest, i, A, edgeSizeThreshold);
    gain += cluster_gain;

    // LOG << "Louvain [" << V(_ipName) << V(_tag) << "]: node " << i << " -> " << A << " (gain: " << gain << ")";
    if (gain > best_gain) {
      best_gain = gain;
      best_partition = A;
    } else if (gain != gain) { // NaN
      if (_context.partition.verbose_output)
        LOG << "Louvain [" << V(_ipName) << V(_tag) << "]: THE GAIN IS NaN: node " << i << " -> " << A << " (gain: " << gain << ")";
    }
  } 

  if (best_gain > eps) {
    // improving = true;
    /// [debug] if(i % 1000 == 0)
    /// [debug] std::cout << "Louvain: node " << i << " -> " << best_partition << " (gain: " << best_gain << ")" << std::endl;
    H_coarsest.changeNodePart(i, part_i, best_partition);
    if (H_coarsest.partSize(part_i) == 0) {
      k_now--;
    }
    return best_gain;
  } else if (! (best_gain <= 0)) {
    if (_context.partition.verbose_output)
      LOG << "Louvain [" << V(_ipName) << V(_tag) << "]: THE GAIN IS TOO SMALL: node " << i << " -> " << best_partition << " (gain: " << best_gain << ")";
  }
  return 0.0;
}

template<typename TypeTraits>
AONCoefficient AONHypermodularityInitialPartitioner<TypeTraits>::QAONGain(
                      SequentialIPHypergraph& H_coarsest, 
                      const HypernodeID i, const PartitionID A, 
                      const HypernodeID edgeSizeThreshold) {
  // Calculate the gain of moving node i to partition A
  // using the AllOrNothing-Hypermodularity-Louvain-Like gain function
  PartitionID part_i = H_coarsest.partID(i);
    ASSERT(A < H_coarsest.k() && A != kInvalidPartition,
         "AONHypermodularityInitialPartitioner::QAONgain: "
         "partition ID " << A << " is invalid");
    ASSERT(part_i < H_coarsest.k() && part_i != kInvalidPartition,
         "AONHypermodularityInitialPartitioner::QAONgain: "
         "partition ID " << part_i  << " is invalid");
  if (part_i == A) {
    return 0.0; // no gain if already in partition A
  }

  AONCoefficient v_A = static_cast<AONCoefficient>(H_coarsest.partOriginalVolume(A));
  AONCoefficient v_i = static_cast<AONCoefficient>(H_coarsest.partOriginalVolume(part_i));
  AONCoefficient d_i = static_cast<AONCoefficient>(H_coarsest.nodeOriginalWeightedDegree(i));

  AONCoefficient delta_vol = 0.0;
  HypernodeID k_max = std::min({H_coarsest.originalMaxEdgeSize(),
                               static_cast<HypernodeID>(H_coarsest.gammaVectorSize() - 1), // zeros at the end of gamma and beta are removed
                               edgeSizeThreshold});
  for (HypernodeID k = 2; k <= k_max; k ++) {
    // _gamma[k] = \beta_k \cdot \gamma_k
    delta_vol += H_coarsest.gamma(k) * (std::pow(v_i, k) - std::pow(v_i - d_i, k) + 
                             std::pow(v_A, k) - std::pow(v_A + d_i, k));
  }

  AONCoefficient delta_cut = 0.0;
  for (const HyperedgeID &he : H_coarsest.incidentEdges(i)) {
    delta_cut += deltaCut(H_coarsest, he, i, part_i, A);
  }
  return delta_cut + delta_vol; // cluster penalty is added outside
}

template<typename TypeTraits>
AONCoefficient AONHypermodularityInitialPartitioner<TypeTraits>::deltaCut(
        const SequentialIPHypergraph& H_coarsest, 
        const HyperedgeID he, const HypernodeID i, const PartitionID part_i, const PartitionID A) {
  // stats needed to distinguish if he is / will be a cutting edge
  HypernodeID size = H_coarsest.edgeSize(he);
  HypernodeID s_he = H_coarsest.originalEdgeSize(he);
  AONCoefficient _beta_S_he = H_coarsest.beta(s_he);
  AONCoefficient weight_he = static_cast<AONCoefficient>(H_coarsest.edgeWeight(he));

  if (size == 1) {
    return 0.0; // single pin nets are never cutting
  }
  
  PartitionID other_node_part = kInvalidPartition;
  bool restNotSame = false;
  for (const HypernodeID &hn : H_coarsest.pins(he)) {
    if (hn != i) {
      if (other_node_part == kInvalidPartition) {
        other_node_part = H_coarsest.partID(hn);
      } else if (other_node_part != H_coarsest.partID(hn)) {
        restNotSame = true;
        break;
      }
    }
  }

  if (restNotSame) {                    
    // other pins are in different parts
    // the edge was and will be cutting
    return 0.0;
  } 
  if (other_node_part == part_i) {
    // all other pins are in part_i
    // the edge was not cutting, but will be after i -> A
    return - _beta_S_he   // _beta[k] = \beta_k
             * weight_he; // as some edges are combined in contraction
  } else if (other_node_part == A) {
    // all other pins are in A
    // the edge was cutting, but won't be after i -> A
    return   _beta_S_he  // _beta[k] = \beta_k
           * weight_he;  // as some edges are combined in contraction
  } else {
    // all other pins are in a third part
    // the edge was and will be cutting
    return 0.0;
  }
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::expandAndCollapse(SequentialIPHypergraph& H_coarsest, std::vector<PartitionID>& z) {
  // Update current communities on H is done in contract(z)
  
  H_coarsest.contract(z /* global community mapping */);
  ASSERT(H_coarsest.isOriginalSizeUsageInParallelNetsDetectionEnabled(), "Original size usage in parallel nets detection is not enabled");
  // return z_changed /* = true */;
}


template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::randomPartition(PartitionedHypergraph& hg) {
  std::uniform_int_distribution<PartitionID> select_random_block(0, _context.partition.k - 1);

  _ip_data.preassignFixedVertices(hg);
  for ( const HypernodeID& hn : hg.nodes() ) {
    if ( !hg.isFixed(hn) ) {
      // Randomly select a block to assign the hypernode
      PartitionID block = select_random_block(_rng);
      PartitionID current_block = block;
      while ( !fitsIntoBlock(hg, hn, current_block) ) {
        // If the hypernode does not fit into the random selected block
        // (because it would violate the balance constraint), we try to
        // assign it to the next block.
        current_block = ( current_block + 1 ) % _context.partition.k;
        if ( current_block == block ) {
          // In case, we find no valid block to assign the current hypernode
          // to, we assign it to random selected block
          break;
        }
      }
      hg.setNodePart(hn, current_block);
    }
  }
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(AONHypermodularityInitialPartitioner)
 
} // namespace mt_kahypar