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

#include "mt-kahypar/definitions.h"
#include "mt-kahypar/utils/randomize.h"
 
/// [debug] #include "mt-kahypar/IFP/mt-kahypar/io/hypergraph_io.h"
 
namespace mt_kahypar {

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::partitionImpl(
    const HypernodeID edgeSizeThreshold, 
    const long long maxNumIter, 
    const /* long */ double eps, 
    const /* long */ double clusterpPenalty,
    const bool randomize,
    bool useOriginalEdgeSizes,
    const InitialPartitioningAlgorithm ipName
    ) {
  /// [debug] std::cout << "partitionImpl()" << std::endl;
  // if num. nodes = k, assign each node to its own block
  // otherwise produce same result as random IP (maybe change later)
  if ( _ip_data.should_initial_partitioner_run(ipName) ) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();

    if (hg.hasFixedVertices() || !hg.hasAON() || !hg.is_static_hypergraph) {
      // fixed vertices are currently not supported in AON-Hypermodularity IP
      randomPartition(hg);
    } else {
      AON_HMLL(hg, 
               edgeSizeThreshold,
               maxNumIter,
               eps,
               clusterpPenalty,
               randomize,
               useOriginalEdgeSizes);
    }
    // =============== General final steps of IP =============

    hg.needsConductancePriorityQueue();

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(ipName, _rng, _tag, time);
  }
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::AON_HMLL(PartitionedHypergraph& hg,
                                                                const HypernodeID edgeSizeThreshold, 
                                                                const long long maxNumIter, 
                                                                const /* long */ double eps, 
                                                                const /* long */ double clusterpPenalty,
                                                                const bool randomize,
                                                                bool useOriginalEdgeSizes) {
  const vec</* long */ double> beta = hg.betaVector();
  const vec</* long */ double> gamma = hg.gammaVector();

  // coarsest underlying hypergraph
  UnderlyingHypergraph H_new = hg.hypergraph().copy();

  if (! useOriginalEdgeSizes || ! H_new.hasOriginalEdgeSizes()) {
    useOriginalEdgeSizes = false;
    if ( ! H_new.hasOriginalEdgeSizes() && useOriginalEdgeSizes)
      if (_context.partition.verbose_output)
        LOG << "No snapshot of original edge size found: using edge sizes of the current hypergraph.";
    // save current edge sizes, weighted degrees and total volume
    H_new.snapshotOriginalEdgeSizes();
    H_new.snapshotOriginalWeightedDegreesAndTotalVolume();
  }
  if (_context.partition.verbose_output)
    if (useOriginalEdgeSizes)
      LOG << "AON IP: Using original edge sizes: # he " << H_new.initialNumEdges() << "# hn " << H_new.initialNumNodes();
    else
      LOG << "AON IP: Using current edge sizes. # he " << H_new.initialNumEdges() << "# hn " << H_new.initialNumNodes();
  H_new.enableSinglePinNetsRemoval(); // single pin nets are never cutting
  H_new.useOriginalSizeInParallelNetsDetection(true); // otherwise gain is incorrect

  // current communities of hg: z: node -> community
  vec<HypernodeID> z(H_new.initialNumNodes(), kInvalidPartition);

  // =====================================================
  //          1. Singleton initial partitioning
  // =====================================================
  for (const HypernodeID &hn : H_new.nodes()) {
    z[hn] = hn;
    H_new.setCommunityID(hn, hn);
  }

  // =====================================================
  //          2. AllOrNothingHMLL: Louvain Cycle
  // =====================================================
  PartitionedHypergraph H_new_partitioned;
  ds::Array<HypernodeID> clusterSizes(H_new.initialNumNodes(), 1);
  vec<HypernodeID> map_z(H_new.initialNumNodes(), kInvalidPartition);
  bool z_changed = false;
  /* long */ double total_gain = 0.0;
  long long counter = 0;
  do {
    counter++;
    /// [debug] std::cout << "Outer Iteration: " << counter << std::endl;
    /** -------------------- Collapse: --------------------
     *  - The community structure on H_new is collapsed by 
     *    merging nodes within the same community;
     *  - H_new_partitioned is rewrited with a singleton
     *    partition on H_new;
     *  - map_z stores the mapping from communityIDs in z
     *    to the HypernodeIDs in H_new which are used as
     *    PartitionIDs in H_new_partitioned.
     */
    collapse(H_new, H_new_partitioned, map_z, clusterSizes);
    /// [debug] std::cout << "Outer Iteration: collapse(..) finished " << counter << std::endl;

    /** ------------------ Louvain Step: ------------------
     *  - Nodes are moved to neighboring partitions as
     *    long as it improves the modularity gain;
     *  - map_z is updated accordingly.
     */
    /* long */ double new_gain = louvainStep(H_new, H_new_partitioned, map_z, clusterSizes,
                                             beta, gamma, edgeSizeThreshold, 
                                             maxNumIter, eps, clusterPenalty, randomize);
    if (_context.partition.verbose_output)
      LOG << "AON IP: step " << counter << " gain = " << new_gain;
    total_gain += new_gain;
    z_changed = (new_gain > eps);
    /// [debug] std::cout << "Outer Iteration: louvainStep(..) finished " << counter << std::endl;

    /** --------------------- Expand: ---------------------
     *  - If H_new_partitioned is still in a singleton 
     *    partition, false is returned;
     *  - otherwise, the community structure on H_new is
     *    updated according to the new partition on 
     *    H_new_partitioned;
     *  - z (communities for H) is updated via map_z;
     *  - true is returned.
     */
    expand(hg, H_new, H_new_partitioned, map_z, z);
    /// [debug] std::cout << "Outer Iteration: expand(..) finished " << counter << std::endl;
  } while (z_changed);

  if (_context.partition.verbose_output)
    LOG << "AON IP finished Louvain with total gain " << total_gain 
        << " and " << H_new_partitioned.k() << " clusters.";

  // =====================================================
  //             3. Finalize Partitioning
  // =====================================================

  // HypernodeID new_k = H_new_partitioned.k();
  // hg.setK(new_k); 
  for (const HypernodeID &hn : hg.nodes()) {
    PartitionID partition = z[hn];
    ASSERT(partition != kInvalidPartition,
        "AONHypermodularityInitialPartitioner::partitionImpl: "
        "node " << hn << " is not assigned to a new partition");
    ASSERT(partition < H_new_partitioned.k(), "AONHypermodularityInitialPartitioner::partitionImpl: "
        "node " << hn << " is assigned to an invalid partition: k = " << H_new_partitioned.k());
    hg.setOnlyNodePart(hn, partition);
  }
  hg.initializePartition();
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::collapse(
        UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, 
        vec<HypernodeID>& map_z, ds::Array<HypernodeID>& clusterSizes) {
  H_new_partitioned = PartitionedHypergraph(H_new.initialNumNodes(), H_new);
  H_new_partitioned.setNecessityOfConductancePriorityQueue(false);
  // create a mapping from the old community IDs to the new ones:
  //     map_z: z-community -> H_new_partitioned-partition
  // (contraction saves old communityID in contracted communityID)
  for (const HypernodeID &hn : H_new.nodes()) {
    H_new_partitioned.setOnlyNodePart(hn, hn);
    clusterSizes[hn] = 1;
    map_z[H_new.communityID(hn)] = hn; // hn is the new partition ID of the node hn
  }
  H_new_partitioned.initializePartition();
}

template<typename TypeTraits>
/* long */ double AONHypermodularityInitialPartitioner<TypeTraits>::louvainStep(
            UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, ds::Array<HypernodeID>& clusterSizes,
            const vec</* long */ double>& beta, const vec</* long */ double>& gamma, const HypernodeID edgeSizeThreshold, 
            const long long maxNumIter, const /* long */ double eps, const /* long */ double clusterpPenalty, const bool randomize) {
  // precompute neighboring nodes
  HypernodeID numNodes = H_new.initialNumNodes();
  ASSERT(H_new_partitioned.k() == static_cast<PartitionID>(numNodes), 
        "AONHypermodularityInitialPartitioner::louvainStep: "
        "no singleton partition at the start: "
        "k = " << H_new_partitioned.k() << ", numNodes = " << numNodes);
  vec<vec<HypernodeID>> neighbors(numNodes, vec<HypernodeID>(0));
  ds::Array<bool> visited(numNodes, false);
  for (const HypernodeID &i : H_new_partitioned.nodes()) {
    visited.assign(numNodes, false);
    for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
      for (const HypernodeID &hn : H_new_partitioned.pins(he)) {
        ASSERT(H_new_partitioned.partID(hn) == static_cast<PartitionID>(hn), 
              "AONHypermodularityInitialPartitioner::louvainStep: "
              "node " << hn << " is not assigned to its singleton partition: "
              "partID(hn) = " << H_new_partitioned.partID(hn));
        visited[hn] = true;
      }
    }
    visited[i] = false;
    for (const HypernodeID& hn : H_new_partitioned.nodes()) {
      if (visited[hn])
        neighbors[i].push_back(hn);
    }

    ASSERT(static_cast<PartitionID>(neighbors[i].size()) < H_new_partitioned.k(), 
           "AONHypermodularityInitialPartitioner::louvainStep: "
           "node " << i << " has more neighboring nodes than partitions");
  }

  /* long */ double total_gain = 0.0;
  bool improving = true;
  long long iter = 0;
  PartitionID k_now = H_new_partitioned.k(); // tracks current number of clusters
  if (randomize) {
    vec<HypernodeID> nodes(numNodes, 0);
    for (HypernodeID i = 0; i < numNodes; ++i) {
      nodes[i] = i;
    }
    std::shuffle(nodes.begin(), nodes.end(), _rng);
    while (improving && (iter++ < maxNumIter)) {
      /* long */ double gain = 0.0;
      improving = false;
      for (const HypernodeID &i : nodes) {
        if ( ! H_new_partitioned.nodeIsEnabled(i) ) continue;
        gain += louvainStepForANode(i, neighbors[i], visited, 
                                       k_now /* reference! */, clusterSizes,
                                       H_new, H_new_partitioned, map_z, 
                                       beta, gamma, edgeSizeThreshold, 
                                       eps, clusterpPenalty, randomize);
      }
      total_gain += gain;
      improving = (gain > eps);
    }
  } else {
    while (improving && (iter++ < maxNumIter)) {
      /* long */ double gain = 0.0;
      improving = false;
      /// [debug]   std::cout << "Louvain: round " << iter << std::endl;
      for (const HypernodeID &i : H_new_partitioned.nodes()) {
        gain += louvainStepForANode(i, neighbors[i], visited, 
                                    k_now /* reference! */, clusterSizes,
                                    H_new, H_new_partitioned, map_z, 
                                    beta, gamma, edgeSizeThreshold, 
                                    eps, clusterpPenalty, randomize);
      }
      total_gain += gain;
      improving = (gain > eps);
    }
  }
  return total_gain;
}

template<typename TypeTraits>
/* long */ double AONHypermodularityInitialPartitioner<TypeTraits>::louvainStepForANode(
        const HypernodeID& i, const vec<HypernodeID>& neighbors_i, ds::Array<bool>& visitedParts, 
        PartitionID& k_now, ds::Array<HypernodeID>& clusterSizes,
        UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, 
        const vec</* long */ double>& beta, const vec</* long */ double>& gamma, const HypernodeID edgeSizeThreshold, 
        const /* long */ double eps, const /* long */ double clusterpPenalty, const bool randomize) {
  unused(randomize);
  /// [debug] if (i % 1000 == 0)
  /// [debug] std::cout << "Louvain: node " << i << std::endl;
  PartitionID part_i = H_new_partitioned.partID(i);

  /// Check all neighboring partitions to find the best gain
  visitedParts.assign(H_new_partitioned.k(), false);
  visitedParts[part_i] = true; // mark current partition as visited

  /* long */ double cluster_gain = 0.0;
  if (clusterpPenalty > 0.0 && clusterSizes[part_i] == 1) {
    cluster_gain = clusterpPenalty * (std::log(k_now) - std::log(k_now - 1));
  }

  /* long */ double best_gain = 0.0;
  PartitionID best_partition = part_i;
  // for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
  //   for (const PartitionID &A : H_new_partitioned.connectivitySet(he)) {
  // bool improving = false;
  for (const HypernodeID &neighbor : neighbors_i) {
    PartitionID A = H_new_partitioned.partID(neighbor);
    if (visitedParts[A]) continue;
    visitedParts[A] = true;
    /* long */ double gain = QAONGain(H_new_partitioned, i, A, beta, gamma, edgeSizeThreshold);
    gain += cluster_gain;

    // LOG << "Louvain: node " << i << " -> " << A << " (gain: " << gain << ")";
    if (gain > best_gain) {
      best_gain = gain;
      best_partition = A;
    } else if (gain != gain) { // NaN
      if (_context.partition.verbose_output)
        LOG << "Louvain: THE GAIN IS NaN: node " << i << " -> " << A << " (gain: " << gain << ")";
    }
  } 

  if (best_gain > eps) {
    // improving = true;
    /// [debug] if(i % 1000 == 0)
    /// [debug] std::cout << "Louvain: node " << i << " -> " << best_partition << " (gain: " << best_gain << ")" << std::endl;
    // Update map_z with the new partition
    ASSERT(part_i == static_cast<PartitionID>(map_z[H_new.communityID(i)]),
             "AONHypermodularityInitialPartitioner::louvainStepForANode: "
             "node " << i << " is not assigned to its mapped community: "
             "partID(i) = " << part_i << ", map_z[communityID(i)] = " << map_z[H_new.communityID(i)]);
    if (H_new_partitioned.changeNodePart(i, part_i, best_partition)) { 
      map_z[H_new.communityID(i)] = best_partition;
      clusterSizes[part_i]--;
      clusterSizes[best_partition]++;
      if (clusterSizes[part_i] == 0) {
        k_now--;
      }
      return best_gain;
    } else {
      LOG << RED << "Louvain: node " << i << " -> " << best_partition << " FAILED";
    }
  } else if (! (best_gain <= 0)) {
    if (_context.partition.verbose_output)
      LOG << "Louvain: THE GAIN IS TOO SMALL: node " << i << " -> " << best_partition << " (gain: " << best_gain << ")";
  }
  return 0.0;
}

template<typename TypeTraits>
/* long */ double AONHypermodularityInitialPartitioner<TypeTraits>::QAONGain(
                      PartitionedHypergraph& H_new_partitioned, 
                      const HypernodeID i, const PartitionID A, 
                      const vec</* long */ double>& beta, const vec</* long */ double>& gamma, const HypernodeID edgeSizeThreshold) {
  unused(beta);
  // Calculate the gain of moving node i to partition A
  // using the AllOrNothing-Hypermodularity-Louvain-Like gain function
  PartitionID part_i = H_new_partitioned.partID(i);
    ASSERT(A < H_new_partitioned.k() && A != kInvalidPartition,
         "AONHypermodularityInitialPartitioner::QAONgain: "
         "partition ID " << A << " is invalid");
    ASSERT(part_i < H_new_partitioned.k() && part_i != kInvalidPartition,
         "AONHypermodularityInitialPartitioner::QAONgain: "
         "partition ID " << part_i  << " is invalid");
  if (part_i == A) {
    return 0.0; // no gain if already in partition A
  }

  /* long */ double v_A = static_cast</* long */ double>(H_new_partitioned.partOriginalVolume(A));
  /* long */ double v_i = static_cast</* long */ double>(H_new_partitioned.partOriginalVolume(part_i));
  /* long */ double d_i = static_cast</* long */ double>(H_new_partitioned.nodeOriginalWeightedDegree(i));

  /* long */ double delta_vol = 0.0;
  HypernodeID k_max = std::min({H_new_partitioned.originalMaxEdgeSize(),
                               static_cast<HypernodeID>(gamma.size() - 1), // zeros at the end of gamma and beta are removed
                               edgeSizeThreshold});
  for (HypernodeID k = 2; k <= k_max; k ++) {
    // _gamma[k] = \beta_k \cdot \gamma_k
    delta_vol += gamma[k] * (std::pow(v_i, k) - std::pow(v_i - d_i, k) + 
                             std::pow(v_A, k) - std::pow(v_A + d_i, k));
  }

  /* long */ double delta_cut = 0.0;
  for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
    // stats needed to distinguish if he is / will be a cutting edge
    HypernodeID pin_count_A = H_new_partitioned.pinCountInPart(he, A);
    HypernodeID pin_count_part_i = H_new_partitioned.pinCountInPart(he, part_i);
    HypernodeID size = H_new_partitioned.edgeSize(he);
    // values needed for the gain computation
    HypernodeID s_he = H_new_partitioned.originalEdgeSize(he);
    /* long */ double _beta_S_he = H_new_partitioned.beta(s_he);
    /* long */ double weight_he = static_cast</* long */ double>(H_new_partitioned.edgeWeight(he));

    // z_he
    if (pin_count_part_i == size && size > 1) {
      // not a cutting edge <=> kroneker_delta(z_he) = 1
      ASSERT(H_new_partitioned.connectivity(he) == 1, 
             "Pin count isn't consistent with connectivity of hyperedge " << he
              << ": connectivity: " << H_new_partitioned.connectivity(he)
              << ": pin count in " << part_i << " : " << pin_count_part_i);
      delta_cut -= _beta_S_he  // _beta[k] = \beta_k
                  * weight_he; // as some edges are combined in contraction
    }
    // z_he i -> A
    if (pin_count_A + 1 == size && size > 1) {
      // won't be a cutting edge after i -> A  <=> kroneker_delta(z_he_i->A) = 1
      ASSERT(H_new_partitioned.connectivity(he) > 1 || size == 1, 
             "Pin count isn't consistent with connectivity of hyperedge " << he
              << ": connectivity: " << H_new_partitioned.connectivity(he)
              << ": pin count in " << A << " : " << pin_count_A);
      delta_cut += _beta_S_he  // _beta[k] = \beta_k
                  * weight_he; // as some edges are combined in contraction
    }
  }
  return delta_cut + delta_vol; // cluster penalty is added outside
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::expand(const PartitionedHypergraph& initPhg, UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, vec<HypernodeID>& z) {

  // Update communities in H_new (to be able to contract it later)
  for (const HypernodeID &hn : H_new.nodes()) {
    HypernodeID partition = H_new_partitioned.partID(hn);
    H_new.setCommunityID(hn, partition);
  }

  // Update current communities on H
  vec<HypernodeID> z_new = z;
  for (const HypernodeID& hn : initPhg.nodes()) {
    HypernodeID community = z[hn];
    ASSERT(static_cast<PartitionID>(map_z[community]) != kInvalidPartition,
            "AONHypermodularityInitialPartitioner::expand : "
            "community " << community << " is not mapped to a new partition");
    z_new[hn] = map_z[community];
  }
  z = z_new;

  vec<HypernodeID> community(H_new.initialNumNodes(), kInvalidPartition);
  for (const HypernodeID& hn : H_new.nodes()) {
    community[hn] = H_new.communityID(hn);
  }
  H_new = H_new.contract(community /* community mapping */);
  ASSERT(H_new.isOriginalSizeUsageInParallelNetsDetectionEnabled(), "Original size usage in parallel nets detection is not enabled");
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