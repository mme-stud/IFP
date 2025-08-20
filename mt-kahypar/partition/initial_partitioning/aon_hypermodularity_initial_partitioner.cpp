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
 
 
namespace mt_kahypar {

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::partitionImpl() {
  // if num. nodes = k, assign each node to its own block
  // otherwise produce same result as random IP (maybe change later)
  if ( _ip_data.should_initial_partitioner_run(InitialPartitioningAlgorithm::aon_hypermodularity) ) {
    HighResClockTimepoint start = std::chrono::high_resolution_clock::now();
    PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();

    if (hg.hasFixedVertices() || !hg.hasAON() || !hg.is_static_hypergraph) {
        // fixed vertices are currently not supported in AON-Hypermodularity IP
        randomPartition(hg);
    } else {
      const vec<double> beta = hg.betaVector();
      const vec<double> gamma = hg.gammaVector();

      // coarsest underlying hypergraph
      UnderlyingHypergraph H = hg.hypergraph().copy();
      
      // save current edge sizes, weighted degrees and total volume
      H.snapshotOriginalEdgeSizes();
      H.snapshotOriginalWeightedDegreesAndTotalVolume();
      H.useOriginalSizeInParallelNetsDetection(true); // otherwise gain is incorrect
      H.disableSinglePinNetsRemoval(); // to not lose contracted edges

      // current communities of hg: z: node -> community
      vec<HypernodeID> z(H.initialNumEdges(), kInvalidPartition);

      // =====================================================
      //          1. Singleton initial partitioning
      // =====================================================
      for (const HypernodeID &hn : H.nodes()) {
        z[hn] = hn;
        H.setCommunityID(hn, hn);
      }

      // =====================================================
      //          2. AllOrNothingHMLL: Louvain Cycle
      // =====================================================
      UnderlyingHypergraph H_new = H.copy();
      PartitionedHypergraph H_new_partitioned;
      vec<HypernodeID> map_z(H.initialNumNodes(), kInvalidPartition);
      bool z_changed = false;
      do {
        /** -------------------- Collapse: --------------------
         *  - The community structure on H_new is collapsed by 
         *    merging nodes within the same community;
         *  - H_new_partitioned is rewrited with a singleton
         *    partition on H_new;
         *  - map_z stores the mapping from communityIDs in z
         *    to the HypernodeIDs in H_new which are used as
         *    PartitionIDs in H_new_partitioned.
         */
        collapse(H_new, H_new_partitioned, map_z);

        /** ------------------ Louvain Step: ------------------
         *  - Nodes are moved to neighboring partitions as
         *    long as it improves the modularity gain;
         *  - map_z is updated accordingly.
         */
        louvainStep(H_new, H_new_partitioned, map_z, beta, gamma);

        /** --------------------- Expand: ---------------------
         *  - If H_new_partitioned is still in a singleton 
         *    partition, false is returned;
         *  - otherwise, the community structure on H_new is
         *    updated according to the new partition on 
         *    H_new_partitioned;
         *  - z (communities for H) is updated via map_z;
         *  - true is returned.
         */
        z_changed = expand(H, H_new, H_new_partitioned, map_z, z);
      } while (z_changed);

      // =====================================================
      //             3. Finalize Partitioning
      // =====================================================

//      HypernodeID new_k = H_new_partitioned.k();
//      hg.setK(new_k, H.initialNumEdges());
      for (const HypernodeID &hn : hg.nodes()) {
        PartitionID partition = z[hn];
        ASSERT(partition != kInvalidPartition,
            "AONHypermodularityInitialPartitioner::partitionImpl: "
            "node " << hn << " is not assigned to a new partition");
        ASSERT(partition < hg.k(), "AONHypermodularityInitialPartitioner::partitionImpl: "
            "node " << hn << " is assigned to an invalid partition: k = " << hg.k());
        hg.setOnlyNodePart(hn, partition);
      }
      hg.initializePartition();
    }

    // =============== General final steps of IP =============

    hg.needsConductancePriorityQueue();

    HighResClockTimepoint end = std::chrono::high_resolution_clock::now();
    double time = std::chrono::duration<double>(end - start).count();
    _ip_data.commit(InitialPartitioningAlgorithm::aon_hypermodularity, _rng, _tag, time);
  }
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::collapse(UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z) {
  vec<HypernodeID> community(H_new.initialNumNodes(), kInvalidPartition);
  for (const HypernodeID& hn : H_new.nodes()) {
    community[hn] = H_new.communityID(hn);
  }
  H_new = H_new.contract(community /* community mapping */);
  H_new_partitioned = PartitionedHypergraph(H_new.initialNumNodes(), H_new);
  H_new_partitioned.setNecessityOfConductancePriorityQueue(false);
  // create a mapping from the old community IDs to the new ones:
  //     map_z: z-community -> H_new_partitioned-partition
  // (contraction saves old communityID in contracted communityID)
  for (const HypernodeID &hn : H_new.nodes()) {
    H_new_partitioned.setOnlyNodePart(hn, hn);
    map_z[H_new.communityID(hn)] = hn; // hn is the new partition ID of the node hn
  }
  H_new_partitioned.initializePartition();
}

template<typename TypeTraits>
void AONHypermodularityInitialPartitioner<TypeTraits>::louvainStep(UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, const vec<double>& beta, const vec<double>& gamma, const long long maxNumIter, const double eps, const bool randomize) {
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

  bool improving = true;
  long long iter = 0;
  while (improving && (iter++ < maxNumIter)) {
    improving = false;
    
    if (randomize) {
      vec<HypernodeID> nodes(numNodes, 0);
      for (HypernodeID i = 0; i < numNodes; ++i) {
        nodes[i] = i;
      }
      std::shuffle(nodes.begin(), nodes.end(), _rng);
      for (const HypernodeID &i : nodes) {
        improving = louvainStepForANode(i, neighbors[i], visited, H_new_partitioned, map_z, beta, gamma, maxNumIter, eps, randomize);
      }
    } else {
      for (const HypernodeID &i : H_new_partitioned.nodes()) {
        improving = louvainStepForANode(i, neighbors[i], visited, H_new_partitioned, map_z, beta, gamma, maxNumIter, eps, randomize);
      }
    }
  }
}

template<typename TypeTraits>
bool AONHypermodularityInitialPartitioner<TypeTraits>::louvainStepForANode(const HypernodeID& i, const vec<HypernodeID>& neighbors_i, const ds::Array<bool>& visitedParts, UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, const vec<double>& beta, const vec<double>& gamma, const long long maxNumIter, const double eps, const bool randomize) {
  HypernodeID part_i = H_new_partitioned.partID(i);

  /// Check all neighboring partitions to find the best gain
  visitedParts.assign(numNodes, false);
  visitedParts[part_i] = true; // mark current partition as visited
  double best_gain = 0.0;
  PartitionID best_partition = part_i;
  // for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
  //   for (const PartitionID &A : H_new_partitioned.connectivitySet(he)) {
  bool improving = false;
  for (const HypernodeID &neighbor : neighbors_i) {
    PartitionID A = H_new_partitioned.partID(neighbor);
    if (visitedParts[A]) continue;
    visitedParts[A] = true;
    double gain = QAONGain(H_new_partitioned, i, A, beta, gamma);
    if (gain > best_gain) {
      best_gain = gain;
      best_partition = A;
    }
  } 

  if (best_gain > eps) {
    improving = true;
    // Update map_z with the new partition
    map_z[H_new.communityID(i)] = best_partition;
    H_new_partitioned.changeNodePart(i, part_i, best_partition);
  }
  return improving;
}

template<typename TypeTraits>
double AONHypermodularityInitialPartitioner<TypeTraits>::QAONGain(PartitionedHypergraph& H_new_partitioned, const HypernodeID i, const PartitionID A, const vec<double>& beta, const vec<double>& gamma) {
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

  double v_A = static_cast<double>(H_new_partitioned.partOriginalVolume(A));
  double v_i = static_cast<double>(H_new_partitioned.partOriginalVolume(part_i));
  double d_i = static_cast<double>(H_new_partitioned.nodeOriginalWeightedDegree(i));

  double delta_vol = 0.0;
  for (HypernodeID k = 1; k <= H_new_partitioned.originalMaxEdgeSize(); k ++) {
    // _gamma[k] = \beta_k \cdot \gamma_k
    delta_vol += gamma[k] * (std::pow(v_i, k) - std::pow(v_i - d_i, k) + 
                             std::pow(v_A, k) - std::pow(v_A + d_i, k));
  }

  double delta_cut = 0.0;
  for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
    // stats needed to distinguish if he is / will be a cutting edge
    HypernodeID pin_count_A = H_new_partitioned.pinCountInPart(he, A);
    HypernodeID pin_count_part_i = H_new_partitioned.pinCountInPart(he, part_i);
    HypernodeID size = H_new_partitioned.edgeSize(he);
    // values needed for the gain computation
    double _beta_S_he = beta[H_new_partitioned.originalEdgeSize(he)];
    double weight_he = static_cast<double>(H_new_partitioned.edgeWeight(he));

    // z_he
    if (pin_count_part_i == size) {
      // not a cutting edge <=> kroneker_delta(z_he) = 1
      ASSERT(H_new_partitioned.connectivity(he) == 1, 
             "Pin count isn't consistent with connectivity of hyperedge " << he
              << ": connectivity: " << H_new_partitioned.connectivity(he)
              << ": pin count in " << part_i << " : " << pin_count_part_i);
      delta_cut -= _beta_S_he  // _beta[k] = \beta_k
                  * weight_he; // as some edges are combined in contraction
    }
    // z_he i -> A
    if (pin_count_A + 1 == size) {
      // won't be a cutting edge after i -> A  <=> kroneker_delta(z_he_i->A) = 1
      ASSERT(H_new_partitioned.connectivity(he) > 1 || size == 1, 
             "Pin count isn't consistent with connectivity of hyperedge " << he
              << ": connectivity: " << H_new_partitioned.connectivity(he)
              << ": pin count in " << A << " : " << pin_count_A);
      delta_cut += _beta_S_he  // _beta[k] = \beta_k
                  * weight_he; // as some edges are combined in contraction
    }
  }

  /* the same, but with 2 additional changeNodePart(..) => 2 loops through incidentEdges 
  // simulate a move i -> A to find new cutting edges
  H_new_partitioned.changeNodePart(i, part_i, A);
  for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
    if (H_new_partitioned.connectivity(he) == 1) { 
      // not a cutting edge <=> kroneker_delta(z_he) = 1
      // minus as _beta[k] = - \beta_k
      delta_cut -= _beta[H_new_partitioned.originalEdgeSize(he)] 
                  * H_new_partitioned.edgeWeight(he);
      }
  }
  // revert the move
  H_new_partitioned.changeNodePart(i, A, part_i);
  for (const HyperedgeID &he : H_new_partitioned.incidentEdges(i)) {
    if (H_new_partitioned.connectivity(he) == 1) { 
      // not a cutting edge <=> kroneker_delta(z_he) = 1
      // plus as _beta[k] = - \beta_k
      delta_cut += _beta[H_new_partitioned.originalEdgeSize(he)]
                  * H_new_partitioned.edgeWeight(he);
    }
  }
  */
  
  return delta_cut + delta_vol;
}

template<typename TypeTraits>
bool AONHypermodularityInitialPartitioner<TypeTraits>::expand(UnderlyingHypergraph& H, UnderlyingHypergraph& H_new, PartitionedHypergraph& H_new_partitioned, vec<HypernodeID>& map_z, vec<HypernodeID>& z) {
  // Check if something changed
  bool z_changed = false;
  vec<bool> notEmptyPart(H_new.initialNumNodes(), false);
  for (const HypernodeID &hn : H_new.nodes()) {
    PartitionID partition = H_new_partitioned.partID(hn);
    ASSERT(partition != kInvalidPartition,
            "AONHypermodularityInitialPartitioner::expand : "
            "partition of hypernode " << hn << " is invalid");
    if (notEmptyPart[partition]) {
      // Not a singleton => changed since the collapse
      z_changed = true;
      break;
    }
    notEmptyPart[partition] = true;
  }
  if (! z_changed)
    return z_changed /* = false */;

  // Update communities in H_new (to be able to contract it later)
  for (const HypernodeID &hn : H_new.nodes()) {
    HypernodeID partition = H_new_partitioned.partID(hn);
    H_new.setCommunityID(hn, partition);
  }

  // Update current communities on H
  vec<HypernodeID> z_new = z;
  for (const HypernodeID &hn : H.nodes()) {
    HypernodeID community = z[hn];
    ASSERT(static_cast<PartitionID>(map_z[community]) != kInvalidPartition,
            "AONHypermodularityInitialPartitioner::expand : "
            "community " << community << " is not mapped to a new partition");
    z_new[hn] = map_z[community];
  }
  z = z_new;
  return z_changed /* = true */;
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