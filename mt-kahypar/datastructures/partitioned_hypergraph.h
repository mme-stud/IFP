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

#pragma once

#include <atomic>
#include <type_traits>
#include <mutex>

#include <tbb/parallel_invoke.h>

#include "kahypar-resources/meta/mandatory.h"

#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/connectivity_info.h"
#include "mt-kahypar/datastructures/streaming_vector.h"
#include "mt-kahypar/datastructures/conductance_pq.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/thread_locals.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {

// Forward
class TargetGraph;

namespace ds {

// Forward
template<typename PartitionedHypergraph, bool maintain_connectivity_set>
class DeltaPartitionedHypergraph;

template <typename Hypergraph = Mandatory,
          typename ConnectivityInformation = ConnectivityInfo>
class PartitionedHypergraph {
 private:
  static_assert(!Hypergraph::is_partitioned,  "Only unpartitioned hypergraphs are allowed");

  using NotificationFunc = std::function<void (SynchronizedEdgeUpdate&)>;
  using DeltaFunction = std::function<void (const SynchronizedEdgeUpdate&)>;
  #define NOOP_NOTIFY_FUNC [] (const SynchronizedEdgeUpdate&) { }
  #define NOOP_FUNC [] (const SynchronizedEdgeUpdate&) { }

  // Factory
  using HypergraphFactory = typename Hypergraph::Factory;

  // REVIEW NOTE: Can't we use a lambda in changeNodePart. And write a second function that calls the first with a lambda that does nothing.
  // Then we could guarantee inlining
  // This would also reduce the code/documentation copy-pasta for with or without gain updates

  static constexpr bool enable_heavy_assert = false;

 public:
  static constexpr bool is_static_hypergraph = Hypergraph::is_static_hypergraph;
  static constexpr bool is_graph = Hypergraph::is_graph;
  static constexpr bool is_partitioned = true;
  static constexpr bool supports_connectivity_set = true;
  static constexpr mt_kahypar_partition_type_t TYPE =
    PartitionedHypergraphType<Hypergraph, ConnectivityInformation>::TYPE;

  static constexpr HyperedgeID HIGH_DEGREE_THRESHOLD = ID(100000);

  using Self = PartitionedHypergraph<Hypergraph, ConnectivityInformation>;
  using UnderlyingHypergraph = Hypergraph;
  using HypernodeIterator = typename Hypergraph::HypernodeIterator;
  using HyperedgeIterator = typename Hypergraph::HyperedgeIterator;
  using IncidenceIterator = typename Hypergraph::IncidenceIterator;
  using IncidentNetsIterator = typename Hypergraph::IncidentNetsIterator;
  using ConInfo = ConnectivityInformation;
  using ConnectivitySetIterator = typename ConnectivityInformation::Iterator;
  template<bool maintain_connectivity_set>
  using DeltaPartition = DeltaPartitionedHypergraph<Self, maintain_connectivity_set>;
  using ExtractedBlock = ExtractedHypergraph<Hypergraph>;

  PartitionedHypergraph() = default;

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph) :
    _input_num_nodes(hypergraph.initialNumNodes()),
    _input_num_edges(hypergraph.initialNumEdges()),
    _k(k),
    _hg(&hypergraph),
    _target_graph(nullptr),
    _conductance_pq(),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_volumes(k, CAtomic<HypernodeWeight>(0)),
    _part_cut_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(
        "Refinement", "part_ids", hypergraph.initialNumNodes(), false, false),
    _con_info(hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize()),
    _pin_count_update_ownership(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true, false) {
    /// [debug] std::cerr << "PartitionedHypergraph::PartitionedHypergraph(k, hypergraph)" << std::endl;
    _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition, false);
  }

  explicit PartitionedHypergraph(const PartitionID k,
                                 Hypergraph& hypergraph,
                                 parallel_tag_t) :
    _input_num_nodes(hypergraph.initialNumNodes()),
    _input_num_edges(hypergraph.initialNumEdges()),
    _k(k),
    _hg(&hypergraph),
    _target_graph(nullptr),
    _conductance_pq(),
    _part_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_volumes(k, CAtomic<HypernodeWeight>(0)),
    _part_cut_weights(k, CAtomic<HypernodeWeight>(0)),
    _part_ids(),
    _con_info(),
    _pin_count_update_ownership() {
    /// [debug] std::cerr << "PartitionedHypergraph::PartitionedHypergraph(k, hypergraph, parallel_tag_t)" << std::endl;
    tbb::parallel_invoke([&] {
      _part_ids.resize(
        "Refinement", "vertex_part_info", hypergraph.initialNumNodes());
      _part_ids.assign(hypergraph.initialNumNodes(), kInvalidPartition);
    }, [&] {
      _con_info = ConnectivityInformation(
        hypergraph.initialNumEdges(), k, hypergraph.maxEdgeSize(), parallel_tag_t { });
    }, [&] {
      _pin_count_update_ownership.resize(
        "Refinement", "pin_count_update_ownership", hypergraph.initialNumEdges(), true);
    });
  }

  // REVIEW NOTE why do we delete copy assignment/construction? wouldn't it be useful to make a copy, e.g. for initial partitioning
  PartitionedHypergraph(const PartitionedHypergraph&) = delete;
  PartitionedHypergraph & operator= (const PartitionedHypergraph &) = delete;

  PartitionedHypergraph(PartitionedHypergraph&& other) = default;
  PartitionedHypergraph & operator= (PartitionedHypergraph&& other) = default;

  ~PartitionedHypergraph() {
    freeInternalData();
  }

  void resetData() {
    /// [debug] std::cerr << "PartitionedHypergraph::resetData()" << std::endl;
    tbb::parallel_invoke([&] {
    }, [&] {
      _part_ids.assign(_part_ids.size(), kInvalidPartition);
    }, [&] {
      _con_info.reset();
    }, [&] {
      for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);
    }, [&] {
      for (auto& x : _part_volumes) x.store(0, std::memory_order_relaxed);
    }, [&] {
      for (auto& x : _part_cut_weights) x.store(0, std::memory_order_relaxed);
    }, [&] {
      if (_conductance_pq.initialized()) {
        _conductance_pq.reset();
      }
    });
  }

  // ####################### General Hypergraph Stats ######################

  Hypergraph& hypergraph() {
    ASSERT(_hg);
    return *_hg;
  }

  void setHypergraph(Hypergraph& hypergraph) {
    /// [debug] std::cerr << "PartitionedHypergraph::setHypergraph(hypergraph)" << std::endl;
    _hg = &hypergraph;
  }

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _hg->initialNumNodes();
  }

  // ! Number of nodes of the input hypergraph
  HypernodeID topLevelNumNodes() const {
    return _input_num_nodes;
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _hg->numRemovedHypernodes();
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _hg->initialNumEdges();
  }

  // ! Number of nodes of the input hypergraph
  HyperedgeID topLevelNumEdges() const {
    return _input_num_edges;
  }

  // ! Number of unique edge ids of the input hypergraph
  HyperedgeID topLevelNumUniqueIds() const {
    return _input_num_edges;
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _hg->initialNumPins();
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _hg->initialTotalVertexDegree();
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _hg->totalWeight();
  }

  // ! Total volume of hypergraph
  HypernodeWeight totalVolume() const {
    return _hg->totalVolume();
  }

  // ! Number of blocks this hypergraph is partitioned into
  PartitionID k() const {
    return _k;
  }


  // ####################### Mapping ######################

  void setTargetGraph(const TargetGraph* target_graph) {
    /// [debug] std::cerr << "PartitionedHypergraph::setTargetGraph(target_graph)" << std::endl;
    _target_graph = target_graph;
  }

  bool hasTargetGraph() const {
    return _target_graph != nullptr;
  }

  const TargetGraph* targetGraph() const {
    return _target_graph;
  }

  // ##################### Conductance  ######################

  // ! Returns the conductance of a block
  double_t conductance(const PartitionID p) const {
    /// [debug] std::cerr << "PartitionedHypergraph::conductance(p)" << std::endl;
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeWeight cut_weight = partCutWeight(p);
    const HypernodeWeight volume = partVolume(p);
    if (volume == 0) {
      return -1;
    }
    return static_cast<double_t>(cut_weight) / static_cast<double_t>(volume);
  }

  // ! Enables the conductance priority queue
  void enableConductancePriorityQueue() {
    /// [debug] std::cerr << "PartitionedHypergraph::enableConductancePriorityQueue()" << std::endl;
    _conductance_pq.initialize(*this);
  }

  bool hasConductancePriorityQueue() const {
    /// [debug] std::cerr << "PartitionedHypergraph::hasConductancePriorityQueue()" << std::endl;
    ASSERT(hasConductancePriorityQueue());
    return _conductance_pq.initialized();
  }

  // ! Get the partition with the hightes conductance
  PartitionID topConductancePart() const {
    /// [debug] std::cerr << "PartitionedHypergraph::topConductancePart()" << std::endl;
    ASSERT(hasConductancePriorityQueue());
    return _conductance_pq.topKey();
  }

  // ! Get the partition with the second hightes conductance
  PartitionID secondTopConductancePart() const {
    /// [debug] std::cerr << "PartitionedHypergraph::secondTopConductancePart()" << std::endl;
    ASSERT(hasConductancePriorityQueue());
    return _conductance_pq.secondTopKey();
  }

  // ! Get top 3 partitions with the highest conductance
  vec<PartitionID> topThreeConductanceParts() const {
    /// [debug] std::cerr << "PartitionedHypergraph::topThreeConductanceParts()" << std::endl;
    ASSERT(hasConductancePriorityQueue());
    return _conductance_pq.topThreeKeys();
  }

  // ! Get a pointer to the conductance priority queue
  const ConductancePriorityQueue<Self>* conductancePriorityQueue() const {
    /// [debug] std::cerr << "PartitionedHypergraph::conductancePriorityQueue()" << std::endl;
    ASSERT(hasConductancePriorityQueue());
    return &_conductance_pq;
  }


  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) {
    /// [debug] std::cerr << "PartitionedHypergraph::doParallelForAllNodes(f)" << std::endl;
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    /// [debug] std::cerr << "PartitionedHypergraph::doParallelForAllNodes(f) const" << std::endl;
    _hg->doParallelForAllNodes(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) {
    /// [debug] std::cerr << "PartitionedHypergraph::doParallelForAllEdges(f)" << std::endl;
    static_cast<const PartitionedHypergraph&>(*this).doParallelForAllEdges(f);
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    /// [debug] std::cerr << "PartitionedHypergraph::doParallelForAllEdges(f) const" << std::endl;
    _hg->doParallelForAllEdges(f);
  }

  // ! Returns an iterator over the set of active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    /// [debug] std::cerr << "PartitionedHypergraph::nodes()" << std::endl;
    return _hg->nodes();
  }

  // ! Returns an iterator over the set of active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    /// [debug] std::cerr << "PartitionedHypergraph::edges()" << std::endl;
    return _hg->edges();
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    /// [debug] std::cerr << "PartitionedHypergraph::incidentEdges(u)" << std::endl;
    return _hg->incidentEdges(u);
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u,
                                                    const size_t pos) const {
    /// [debug] std::cerr << "PartitionedHypergraph::incidentEdges(u, pos)" << std::endl;
    return _hg->incident_nets_of(u, pos);
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    /// [debug] std::cerr << "PartitionedHypergraph::pins(e)" << std::endl;
    return _hg->pins(e);
  }

  // ! Returns a range to loop over the set of block ids contained in hyperedge e.
  IteratorRange<ConnectivitySetIterator> connectivitySet(const HyperedgeID e) const {
    /// [debug] std::cerr << "PartitionedHypergraph::connectivitySet(e)" << std::endl;
    ASSERT(_hg->edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    return _con_info.connectivitySet(e);
  }

  // ####################### Hypernode Information #######################

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return _hg->nodeWeight(u);
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    /// [debug] std::cerr << "PartitionedHypergraph::setNodeWeight(u, weight)" << std::endl;
    const PartitionID block = partID(u);
    if ( block != kInvalidPartition ) {
      ASSERT(block < _k);
      const HypernodeWeight delta = weight - _hg->nodeWeight(u);
      _part_weights[block] += delta;
    }
    _hg->setNodeWeight(u, weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    return _hg->nodeDegree(u);
  }

  // ! Weighted degree od a hypernode
  HyperedgeWeight nodeWeightedDegree(const HypernodeID u) const {
    return _hg->nodeWeightedDegree(u);
  }

  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return _hg->nodeIsEnabled(u);
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    /// [debug] std::cerr << "PartitionedHypergraph::enableHypernode(u)" << std::endl;
    _hg->enableHypernode(u);
  }

  // ! Disable a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    /// [debug] std::cerr << "PartitionedHypergraph::disableHypernode(u)" << std::endl;
    _hg->disableHypernode(u);
  }

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u, const PartitionID to) {
    /// [debug] std::cerr << "PartitionedHypergraph::restoreDegreeZeroHypernode(u, to)" << std::endl;
    _hg->restoreDegreeZeroHypernode(u);
    setNodePart(u, to);
  }

  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  HypernodeWeight edgeWeight(const HyperedgeID e) const {
    return _hg->edgeWeight(e);
  }

  // ! Unique id of a hyperedge
  HyperedgeID uniqueEdgeID(const HyperedgeID e) const {
    return e;
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    /// [debug] std::cerr << "PartitionedHypergraph::setEdgeWeight(e, weight)" << std::endl;
    _hg->setEdgeWeight(e, weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    return _hg->edgeSize(e);
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return _hg->edgeIsEnabled(e);
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    /// [debug] std::cerr << "PartitionedHypergraph::enableHyperedge(e)" << std::endl;
    _hg->enableHyperedge(e);
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    /// [debug] std::cerr << "PartitionedHypergraph::disableHyperedge(e)" << std::endl;
    _hg->disableHyperedge(e);
  }

    // ! Target of an edge
  HypernodeID edgeTarget(const HyperedgeID) const {
    throw UnsupportedOperationException(
      "edgeTarget(e) is only supported on graph data structure");
    return kInvalidHypernode;
  }

  // ! Source of an edge
  HypernodeID edgeSource(const HyperedgeID) const {
    throw UnsupportedOperationException(
      "edgeSource(e) is only supported on graph data structure");
    return kInvalidHypernode;
  }

  // ! Whether the edge is a single pin edge
  bool isSinglePin(const HyperedgeID) const {
    throw UnsupportedOperationException(
      "isSinglePin(e) is only supported on graph data structure");
    return false;
  }

  // ####################### Uncontraction #######################

  /**
   * Uncontracts a batch of contractions in parallel. The batches must be uncontracted exactly
   * in the order computed by the function createBatchUncontractionHierarchy(...).
   */
  template<typename GainCache>
  void uncontract(const Batch& batch, GainCache& gain_cache) {
    /// [debug] std::cerr << "PartitionedHypergraph::uncontract(batch, gain_cache)" << std::endl;
    // Set block ids of contraction partners
    tbb::parallel_for(UL(0), batch.size(), [&](const size_t i) {
      const Memento& memento = batch[i];
      ASSERT(nodeIsEnabled(memento.u));
      ASSERT(!nodeIsEnabled(memento.v));
      const PartitionID part_id = partID(memento.u);
      ASSERT(part_id != kInvalidPartition && part_id < _k);
      setOnlyNodePart(memento.v, part_id);
      // update _part_volumes part 1
      decrementVolumeOfBlock(part_id, nodeWeightedDegree(memento.u));
    });

    _hg->uncontract(batch,
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
        // In this case, u and v are incident to hyperedge he after uncontraction
        const PartitionID block = partID(u);
        const HypernodeID pin_count_in_part_after = incrementPinCountOfBlock(he, block);
        ASSERT(pin_count_in_part_after > 1, V(u) << V(v) << V(he));
        gain_cache.uncontractUpdateAfterRestore(*this, u, v, he, pin_count_in_part_after);
      },
      [&](const HypernodeID u, const HypernodeID v, const HyperedgeID he) {
        gain_cache.uncontractUpdateAfterReplacement(*this, u, v, he);
      }
    );
    
    // update _part_volumes part 2
    tbb::parallel_for(UL(0), batch.size(), [&](const size_t i) {
      const Memento& memento = batch[i];
      ASSERT(nodeIsEnabled(memento.u));
      ASSERT(nodeIsEnabled(memento.v));
      const PartitionID part_id = partID(memento.u);
      ASSERT(part_id != kInvalidPartition && part_id < _k && part_id == partID(memento.v));
      incrementVolumeOfBlock(part_id, nodeWeightedDegree(memento.u) + nodeWeightedDegree(memento.v));
    });

    // update _conductance_pq if enabled  
    // after this the gain cache should be updated (?)
    // => ConductanceGainCache::initializes_gain_cache_entry_after_batch_uncontractions = true ?
    if (hasConductancePriorityQueue()) {
      _conductance_pq.globalUpdate(*this);
    }

    if constexpr ( GainCache::initializes_gain_cache_entry_after_batch_uncontractions ) {
      tbb::parallel_for(UL(0), batch.size(), [&](const size_t i) {
        const Memento& memento = batch[i];
        gain_cache.initializeGainCacheEntryForNode(*this, memento.v);
      });
    }
    gain_cache.batchUncontractionsCompleted();
  }

  // ####################### Restore Hyperedges #######################

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    /// [debug] std::cerr << "PartitionedHypergraph::restoreLargeEdge(he)" << std::endl;
    _hg->restoreLargeEdge(he);

    // Recalculate pin count in parts
    // and part volumes
    const size_t incidence_array_start = _hg->hyperedge(he).firstEntry();
    const size_t incidence_array_end = _hg->hyperedge(he).firstInvalidEntry();
    tbb::enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);
    tbb::enumerable_thread_specific< vec<HyperedgeWeight> > ets_add_part_volumes(_k, 0);
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _hg->_incidence_array[pos];
      const PartitionID block = partID(pin);
      ++ets_pin_count_in_part.local()[block];
      ets_add_part_volumes.local()[block] += edgeWeight(he);//  * ets_pin_count_in_part.local()[block]; - wrong: for each pin only once
    });

    // Aggregate local additional _part_volumes for each block
    for ( PartitionID block = 0; block < _k; ++block ) {
      HyperedgeWeight add_part_volume = 0;
      for ( const vec<HyperedgeWeight>& local_add_part_volumes : ets_add_part_volumes) {
        add_part_volume += local_add_part_volumes[block];
      }
      incrementVolumeOfBlock(block, add_part_volume);
    }

    // Aggregate local pin count for each block
    // and update _part_cut_weights
    for ( PartitionID block = 0; block < _k; ++block ) {
      HypernodeID pin_count_in_part = 0;
      for ( const vec<HypernodeID>& local_pin_count : ets_pin_count_in_part ) {
        pin_count_in_part += local_pin_count[block];
      }

      if ( pin_count_in_part > 0 ) {
        _con_info.setPinCountInPart(he, block, pin_count_in_part);
        _con_info.addBlock(he, block);
        // update _part_cut_weights
        if (pin_count_in_part < edgeSize(he)) {
          incrementCutWeightOfBlock(block, edgeWeight(he));
        }
      }
    }

    // update _conductance_pq for changed partitions (if enabled)
    if (hasConductancePriorityQueue()) {
      _conductance_pq.globalUpdate(*this, he);
    }
  }

  /**
   * Restores a previously removed set of singple-pin and parallel hyperedges. Note, that hes_to_restore
   * must be exactly the same and given in the reverse order as returned by removeSinglePinAndParallelNets(...).
   */
  template<typename GainCache>
  void restoreSinglePinAndParallelNets(const vec<typename Hypergraph::ParallelHyperedge>& hes_to_restore,
                                       GainCache& gain_cache) {
    /// [debug] std::cerr << "PartitionedHypergraph::restoreSinglePinAndParallelNets(hes_to_restore, gain_cache)" << std::endl;
    // needed decide if conductance_pq should be updated
    HyperedgeWeight old_total_volume = totalVolume();

    // Restore hyperedges in hypergraph
    _hg->restoreSinglePinAndParallelNets(hes_to_restore);

    // Compute pin counts of restored hyperedges and gain cache values of vertices contained
    // single-pin hyperedges. Note, that restoring parallel hyperedges does not change any
    // value in the gain cache, since it already contributes to the gain via its representative.
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);
    // TODO: ets_pin_count_in_part - unused???
    // also compute additional volumes for all parts: only single-pin nets bring new volume
    //      (_part_volume is atomic => no need for tls_enumerable_thread_specific)
    tbb::parallel_for(UL(0), hes_to_restore.size(), [&](const size_t i) {
      const HyperedgeID he = hes_to_restore[i].removed_hyperedge;
      const HyperedgeID representative = hes_to_restore[i].representative;
      ASSERT(edgeIsEnabled(he));
      const bool is_single_pin_he = edgeSize(he) == 1;
      if ( is_single_pin_he ) {
        // Restore single-pin net
        HypernodeID single_vertex_of_he = kInvalidHypernode;
        for ( const HypernodeID& pin : pins(he) ) {
          single_vertex_of_he = pin;
        }
        ASSERT(single_vertex_of_he != kInvalidHypernode);

        const PartitionID block_of_single_pin = partID(single_vertex_of_he);
        _con_info.addBlock(he, block_of_single_pin);
        _con_info.setPinCountInPart(he, block_of_single_pin, 1);
        gain_cache.restoreSinglePinHyperedge(
          single_vertex_of_he, block_of_single_pin, edgeWeight(he));
        _part_volumes[block_of_single_pin] += edgeWeight(he);
      } else {
        // Restore parallel net => pin count information given by representative
        ASSERT(edgeIsEnabled(representative));
        for ( const PartitionID& block : connectivitySet(representative) ) {
          _con_info.addBlock(he, block);
          _con_info.setPinCountInPart(he, block, pinCountInPart(representative, block));
        }
        gain_cache.restoreIdenticalHyperedge(*this, he);

        HEAVY_REFINEMENT_ASSERT([&] {
          for ( PartitionID block = 0; block < _k; ++block ) {
            if ( pinCountInPart(he, block) != pinCountInPartRecomputed(he, block) ) {
              LOG << "Pin count in part of hyperedge" << he << "in block" << block
                  << "is" << pinCountInPart(he, block) << ", but should be"
                  << pinCountInPartRecomputed(he, block);
              return false;
            }
          }
          return true;
        }());
      }
    });

    // update _conductance_pq (if enabled)
    if (hasConductancePriorityQueue()) {
      /*  Only if the total volume of the hypergraph has changed, 
       * we need to update the conductance priority queue
       * (<=> if at least one single-pin net was restored).
       */
      if (totalVolume() != old_total_volume) {
        _conductance_pq.globalUpdate(*this);
      }
    }
  }

  // ####################### Partition Information #######################

  // ! Block that vertex u belongs to
  PartitionID partID(const HypernodeID u) const {
    ASSERT(u < initialNumNodes(), "Hypernode" << u << "does not exist");
    return _part_ids[u];
  }

  void extractPartIDs(Array<PartitionID>& part_ids) {
    /// [debug] std::cerr << "PartitionedHypergraph::extractPartIDs(part_ids)" << std::endl;
    // If we pass the input hypergraph to initial partitioning, then initial partitioning
    // will pass an part ID vector of size |V'|, where V' are the number of nodes of
    // smallest hypergraph, while the _part_ids vector of the input hypergraph is initialized
    // with the original number of nodes. This can cause segmentation fault when we simply swap them
    // during main uncoarsening.
    if ( _part_ids.size() == part_ids.size() ) {
      std::swap(_part_ids, part_ids);
    } else {
      ASSERT(part_ids.size() <= _part_ids.size());
      tbb::parallel_for(UL(0), part_ids.size(), [&](const size_t i) {
        part_ids[i] = _part_ids[i];
      });
    }
  }

  void setOnlyNodePart(const HypernodeID u, PartitionID p) {
    /// [debug] std::cerr << "PartitionedHypergraph::setOnlyNodePart(u, p)" << std::endl;
    ASSERT(p != kInvalidPartition && p < _k);
    ASSERT(_part_ids[u] == kInvalidPartition);
    _part_ids[u] = p;
  }

  // ! Sets a partition of a node for the first time. 
  // ! This function is surely not thread-safe after adding _part_cut_weights 
  // ! (mostly used for tests?)
  void setNodePart(const HypernodeID u, PartitionID p) {
    /// [debug] std::cerr << "PartitionedHypergraph::setNodePart(u, p)" << std::endl;
    setOnlyNodePart(u, p);
    _part_weights[p].fetch_add(nodeWeight(u), std::memory_order_relaxed);
    incrementVolumeOfBlock(p, nodeWeightedDegree(u));
    for (HyperedgeID he : incidentEdges(u)) {
      const PartitionID old_connectivity = connectivity(he);
      incrementPinCountOfBlock(he, p);
      const PartitionID new_connectivity = connectivity(he);
      // update _part_cut_weights
      if (old_connectivity == 1 && new_connectivity > 1) {
        // Hyperedge he is a new cutting edge for all 2 partitions of its pins
        ASSERT(new_connectivity == 2);
        for (PartitionID partId : connectivitySet(he)) {
          incrementCutWeightOfBlock(partId, edgeWeight(he));
        }
      } else if (new_connectivity != old_connectivity && new_connectivity > 1) {
        ASSERT(new_connectivity == old_connectivity + 1);
        // Hyperedge he was a cutting edge for partitions of its pins
        // and is now a cutting edge for a new partition p of its pin u
        // _part_cut_weights[p].fetch_add(edgeWeight(he), std::memory_order_relaxed);
        incrementCutWeightOfBlock(p, edgeWeight(he));
      }
    }
  }

  // ! Changes the block id of vertex u from block 'from' to block 'to'
  // ! Returns true, if move of vertex u to corresponding block succeeds.
  template<typename SuccessFunc>
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      HypernodeWeight max_weight_to,
                      SuccessFunc&& report_success,
                      const DeltaFunction& delta_func,
                      const NotificationFunc& notify_func = NOOP_NOTIFY_FUNC,
                      const bool force_moving_fixed_vertices = false) {
    /// [debug] std::cerr << "PartitionedHypergraph::changeNodePart(u, from, to, max_weight_to, report_success, delta_func, notify_func, force_moving_fixed_vertices)" << std::endl;
    unused(force_moving_fixed_vertices);
    ASSERT(partID(u) == from);
    ASSERT(from != to);
    ASSERT(force_moving_fixed_vertices || !isFixed(u));
    const HypernodeWeight wu = nodeWeight(u);
    const HypernodeWeight to_weight_after = _part_weights[to].add_fetch(wu, std::memory_order_relaxed);
    if (to_weight_after <= max_weight_to) {
      _part_ids[u] = to;
      _part_weights[from].fetch_sub(wu, std::memory_order_relaxed);
      // update _part_volumes
      decrementVolumeOfBlock(from, nodeWeightedDegree(u));
      incrementVolumeOfBlock(to, nodeWeightedDegree(u));
      report_success();
      SynchronizedEdgeUpdate sync_update;
      sync_update.from = from;
      sync_update.to = to;
      sync_update.target_graph = _target_graph;
      sync_update.edge_locks = &_pin_count_update_ownership;
      for ( const HyperedgeID he : incidentEdges(u) ) {
        // updates _part_cut_weights in updatePinCountOfHyperedge(...)
        updatePinCountOfHyperedge(he, from, to, sync_update, delta_func, notify_func);
        // TODO (?): SPLIT INTO TWO FUNCTIONS? -> no... We need to update _part_cut_weights behind the lock
      }
      // update _conductance_pq if enabled: do it after updating _part_cut_weights and _part_volumes
      if (hasConductancePriorityQueue()) {
        _conductance_pq.adjustKey(from, partCutWeight(from), partVolume(from));
        _conductance_pq.adjustKey(to, partCutWeight(to), partVolume(to));
      }
      return true;
    } else {
      _part_weights[to].fetch_sub(wu, std::memory_order_relaxed);
      return false;
    }
  }

  // curry
  bool changeNodePart(const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      const DeltaFunction& delta_func = NOOP_FUNC,
                      const bool force_moving_fixed_vertex = false) {
    return changeNodePart(u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{},
        delta_func, NOOP_NOTIFY_FUNC, force_moving_fixed_vertex);
  }

  template<typename GainCache, typename SuccessFunc>
  bool changeNodePart(GainCache& gain_cache,
                      const HypernodeID u,
                      PartitionID from,
                      PartitionID to,
                      HypernodeWeight max_weight_to,
                      SuccessFunc&& report_success,
                      const DeltaFunction& delta_func) {
    auto my_delta_func = [&](const SynchronizedEdgeUpdate& sync_update) {
      delta_func(sync_update);
      gain_cache.deltaGainUpdate(*this, sync_update);
    };
    if constexpr ( !GainCache::requires_notification_before_update ) {
      return changeNodePart(u, from, to, max_weight_to, report_success, my_delta_func);
    } else {
      return changeNodePart(u, from, to, max_weight_to, report_success, my_delta_func,
        [&](SynchronizedEdgeUpdate& sync_update) {
          sync_update.pin_count_in_from_part_after = pinCountInPart(sync_update.he, from) - 1;
          sync_update.pin_count_in_to_part_after = pinCountInPart(sync_update.he, to) + 1;
          gain_cache.notifyBeforeDeltaGainUpdate(*this, sync_update);
        });
    }
  }

  template<typename GainCache>
  bool changeNodePart(GainCache& gain_cache,
                      const HypernodeID u,
                      PartitionID from,
                      PartitionID to) {
    return changeNodePart(gain_cache, u, from, to,
      std::numeric_limits<HypernodeWeight>::max(), []{}, NoOpDeltaFunc());
  }

  // ! Weight of a block
  HypernodeWeight partWeight(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_weights[p].load(std::memory_order_relaxed);
  }
  
  // ! Volume of a block
  HyperedgeWeight partVolume(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_volumes[p].load(std::memory_order_relaxed);
  }

  // ! Cut weight of a block 
  HypernodeWeight partCutWeight(const PartitionID p) const {
    ASSERT(p != kInvalidPartition && p < _k);
    return _part_cut_weights[p].load(std::memory_order_relaxed);
  }

  // ! Returns, whether hypernode u is adjacent to a least one cut hyperedge.
  bool isBorderNode(const HypernodeID u) const {
    if ( nodeDegree(u) <= HIGH_DEGREE_THRESHOLD ) {
      for ( const HyperedgeID& he : incidentEdges(u) ) {
        if ( connectivity(he) > 1 ) {
          return true;
        }
      }
      return false;
    } else {
      // TODO maybe we should allow these in label propagation? definitely not in FM
      // In case u is a high degree vertex, we omit the border node check and
      // and return false. Assumption is that it is very unlikely that such a
      // vertex can change its block.
      return false;
    }
  }

  HypernodeID numIncidentCutHyperedges(const HypernodeID u) const {
    HypernodeID num_incident_cut_hyperedges = 0;
    for ( const HyperedgeID& he : incidentEdges(u) ) {
      if ( connectivity(he) > 1 ) {
        ++num_incident_cut_hyperedges;
      }
    }
    return num_incident_cut_hyperedges;
  }

  // ! Number of blocks which pins of hyperedge e belongs to
  PartitionID connectivity(const HyperedgeID e) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    return _con_info.connectivity(e);
  }

  // ! Returns the number pins of hyperedge e that are part of block id
  HypernodeID pinCountInPart(const HyperedgeID e, const PartitionID p) const {
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    return _con_info.pinCountInPart(e, p);
  }

  // ! Creates a shallow copy of the connectivity set of hyperedge he
  StaticBitset& shallowCopyOfConnectivitySet(const HyperedgeID he) const {
    return _con_info.shallowCopy(he);
  }

  // ! Creates a deep copy of the connectivity set of hyperedge he
  Bitset& deepCopyOfConnectivitySet(const HyperedgeID he) const {
    return _con_info.deepCopy(he);
  }

  const ConInfo& getConnectivityInformation() const {
    return _con_info;
  }

  // ! Initializes the partition of the hypergraph, if block ids are assigned with
  // ! setOnlyNodePart(...). In that case, block weights and pin counts in part for
  // ! each hyperedge must be initialized explicitly here.
  void initializePartition() {
    /// [debug] std::cerr << "PartitionedHypergraph::initializePartition()" << std::endl;
    tbb::parallel_invoke(
            [&] { initializeBlockWeights(); },
            [&] { initializeBlockVolumes(); },
            [&] { 
              initializePinCountInPart();
              initializeBlockCutWeights(); // must be called after initializePinCountInPart()
            }
    );
  }

  // ! Reset partition (not thread-safe)
  void resetPartition() {
    /// [debug] std::cerr << "PartitionedHypergraph::resetPartition()" << std::endl;
    _part_ids.assign(_part_ids.size(), kInvalidPartition, false);
    for (auto& x : _part_weights) x.store(0, std::memory_order_relaxed);
    for (auto& x : _part_volumes) x.store(0, std::memory_order_relaxed);
    for (auto& x : _part_cut_weights) x.store(0, std::memory_order_relaxed);
    _conductance_pq.reset();

    // Reset pin count in part and connectivity set
    _con_info.reset(false);
  }

  // ! Only for testing
  void recomputePartWeights() {
    /// [debug] std::cerr << "PartitionedHypergraph::recomputePartWeights()" << std::endl;
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].store(0);
    }

    for (HypernodeID u : nodes()) {
      _part_weights[ partID(u) ] += nodeWeight(u);
    }
  }

  // ! Only for testing
  void recomputePartVolumes() {
    /// [debug] std::cerr << "PartitionedHypergraph::recomputePartVolumes()" << std::endl;
    for (PartitionID p = 0; p < _k; ++p) {
      _part_volumes[p].store(0);
    }

    for (HypernodeID u : nodes()) {
      for (const HyperedgeID& e : incidentEdges(u)) {
        _part_volumes[ partID(u) ] += edgeWeight(e);
      }
    }
  }

  // ! Only for testing
  void recomputePartCutWeights() {
    /// [debug] std::cerr << "PartitionedHypergraph::recomputePartCutWeights()" << std::endl;
    for (PartitionID p = 0; p < _k; ++p) {
      _part_cut_weights[p].store(0);
    }

    for (HyperedgeID he : edges()) {
      if ( connectivity(he) > 1 ) {
        for ( PartitionID block : connectivitySet(he) ) {
          _part_cut_weights[block] += edgeWeight(he);
        }
      }
    }
  }

  // ! Only for testing
  void recomputeConductancePriorityQueue() {
    /// [debug] std::cerr << "PartitionedHypergraph::recomputeConductancePriorityQueue()" << std::endl;
    if (hasConductancePriorityQueue()) {
      _conductance_pq.globalUpdate(*this);
    }
  }

  // ! Only for testing
  bool checkConductancePriorityQueue() {
    /// [debug] std::cerr << "PartitionedHypergraph::checkConductancePriorityQueue()" << std::endl;
    if (hasConductancePriorityQueue()) {
      return _conductance_pq.check(*this);
    }
    return true;
  }

  // ! Only for testing
  bool checkTrackedPartitionInformation() {
    /// [debug] std::cerr << "PartitionedHypergraph::checkTrackedPartitionInformation()" << std::endl;
    bool success = true;

    for (HyperedgeID e : edges()) {
      PartitionID expected_connectivity = 0;
      for (PartitionID i = 0; i < k(); ++i) {
        const HypernodeID actual_pin_count_in_part = pinCountInPart(e, i);
        if ( actual_pin_count_in_part != pinCountInPartRecomputed(e, i) ) {
          LOG << "Pin count of hyperedge" << e << "in block" << i << "=>" <<
              "Expected:" << V(pinCountInPartRecomputed(e, i)) << "," <<
              "Actual:" <<  V(pinCountInPart(e, i));
          success = false;
        }
        expected_connectivity += (actual_pin_count_in_part > 0);
      }
      if ( expected_connectivity != connectivity(e) ) {
        LOG << "Connectivity of hyperedge" << e << "=>" <<
            "Expected:" << V(expected_connectivity)  << "," <<
            "Actual:" << V(connectivity(e));
        success = false;
      }
    }
    return success;
  }

  // ! Only for testing
  template<typename GainCache>
  bool checkTrackedPartitionInformation(GainCache& gain_cache) {
    /// [debug] std::cerr <<  "PartitionedHypergraph::checkTrackedPartitionInformation(gain_cache)" << std::endl;
    bool success = true;

    for (HyperedgeID e : edges()) {
      PartitionID expected_connectivity = 0;
      for (PartitionID i = 0; i < k(); ++i) {
        const HypernodeID actual_pin_count_in_part = pinCountInPart(e, i);
        if ( actual_pin_count_in_part != pinCountInPartRecomputed(e, i) ) {
          LOG << "Pin count of hyperedge" << e << "in block" << i << "=>" <<
              "Expected:" << V(pinCountInPartRecomputed(e, i)) << "," <<
              "Actual:" <<  V(pinCountInPart(e, i));
          success = false;
        }
        expected_connectivity += (actual_pin_count_in_part > 0);
      }
      if ( expected_connectivity != connectivity(e) ) {
        LOG << "Connectivity of hyperedge" << e << "=>" <<
            "Expected:" << V(expected_connectivity)  << "," <<
            "Actual:" << V(connectivity(e));
        success = false;
      }
    }

    if ( gain_cache.isInitialized() ) {
      for (HypernodeID u : nodes()) {
        const PartitionID block_of_u = partID(u);
        if ( gain_cache.penaltyTerm(u, block_of_u) !=
             gain_cache.recomputePenaltyTerm(*this, u) ) {
          LOG << "Penalty term of hypernode" << u << "=>" <<
              "Expected:" << V(gain_cache.recomputePenaltyTerm(*this, u)) << ", " <<
              "Actual:" <<  V(gain_cache.penaltyTerm(u, block_of_u));
          for ( const HyperedgeID& e : incidentEdges(u) ) {
            LOG << V(u) << V(partID(u)) << V(e) << V(edgeSize(e))
                << V(edgeWeight(e)) << V(pinCountInPart(e, partID(u)));
          }
          success = false;
        }

        for (const PartitionID& i : gain_cache.adjacentBlocks(u)) {
          if (partID(u) != i) {
            if ( gain_cache.benefitTerm(u, i) !=
                 gain_cache.recomputeBenefitTerm(*this, u, i) ) {
              LOG << "Benefit term of hypernode" << u << "in block" << i << "=>" <<
                  "Expected:" << V(gain_cache.recomputeBenefitTerm(*this, u, i)) << ", " <<
                  "Actual:" <<  V(gain_cache.benefitTerm(u, i));
              success = false;
            }
          }
        }
      }
      if ( !gain_cache.verifyTrackedAdjacentBlocksOfNodes(*this) ) {
        success = false;
      }
    }
    return success;
  }

  // ####################### Fixed Vertex Support #######################

  bool hasFixedVertices() const {
    return _hg->hasFixedVertices();
  }

  bool isFixed(const HypernodeID hn) const {
    return _hg->isFixed(hn);
  }

  PartitionID fixedVertexBlock(const HypernodeID hn) const {
    return _hg->fixedVertexBlock(hn);
  }

  // ####################### Memory Consumption #######################

  // ! No info about the memory consumption of the conductance priority queue
  void memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);

    utils::MemoryTreeNode* hypergraph_node = parent->addChild("Hypergraph");
    _hg->memoryConsumption(hypergraph_node);
    utils::MemoryTreeNode* connectivity_info_node = parent->addChild("Connectivity Information");
    _con_info.memoryConsumption(connectivity_info_node);

    parent->addChild("Part Weights", sizeof(CAtomic<HypernodeWeight>) * _k);
    parent->addChild("Part IDs", sizeof(PartitionID) * _hg->initialNumNodes());
    parent->addChild("HE Ownership", sizeof(SpinLock) * _hg->initialNumNodes());
  }

  // ####################### Extract Block #######################

  // ! Extracts a block of the partition.
  // ! The extracted block stores a hypergraph containing all nodes and hyperedges
  // ! of the corresponding block, a vector that maps the node IDs of the original
  // ! hypergraph to the node IDs of the extracted hypergraph, and an array that
  // ! stores which hyperedges of the extracted hypergraph are already cut in
  // ! the original hypergraph.
  // ! If cut_net_splitting is activated, then cut hyperedges are splitted containing
  // ! only the pins of the corresponding block. Otherwise, they are discarded.
  ExtractedBlock extract(const PartitionID block,
                         const vec<uint8_t>* already_cut,
                         bool cut_net_splitting,
                         bool stable_construction_of_incident_edges) {
    /// [debug] std::cerr << "PartitionedHypergraph::extract(block, already_cut, cut_net_splitting, stable_construction_of_incident_edges)" << std::endl;
    ASSERT(block != kInvalidPartition && block < _k);
    ASSERT(!already_cut || already_cut->size() == _hg->initialNumEdges());

    // Compactify vertex ids
    ExtractedBlock extracted_block;
    vec<HypernodeID>& hn_mapping = extracted_block.hn_mapping;
    hn_mapping.assign(_hg->initialNumNodes(), kInvalidHypernode);
    vec<HyperedgeID> he_mapping(_hg->initialNumEdges(), kInvalidHyperedge);
    HypernodeID num_hypernodes = 0;
    HypernodeID num_hyperedges = 0;
    tbb::parallel_invoke([&] {
      for ( const HypernodeID& hn : nodes() ) {
        if ( partID(hn) == block ) {
          hn_mapping[hn] = num_hypernodes++;
        }
      }
    }, [&] {
      for ( const HyperedgeID& he : edges() ) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          he_mapping[he] = num_hyperedges++;
        }
      }
    });

    // Extract plain hypergraph data for corresponding block
    using HyperedgeVector = vec<vec<HypernodeID>>;
    HyperedgeVector edge_vector;
    vec<HyperedgeWeight> hyperedge_weight;
    vec<HypernodeWeight> hypernode_weight;
    vec<uint8_t> extracted_already_cut;
    tbb::parallel_invoke([&] {
      edge_vector.resize(num_hyperedges);
      hyperedge_weight.resize(num_hyperedges);
      doParallelForAllEdges([&](const HyperedgeID he) {
        if ( pinCountInPart(he, block) > 1 &&
             (cut_net_splitting || connectivity(he) == 1) ) {
          ASSERT(he_mapping[he] < num_hyperedges);
          hyperedge_weight[he_mapping[he]] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == block ) {
              edge_vector[he_mapping[he]].push_back(hn_mapping[pin]);
            }
          }
        }
      });
    }, [&] {
      hypernode_weight.resize(num_hypernodes);
      doParallelForAllNodes([&](const HypernodeID hn) {
        if ( partID(hn) == block ) {
          hypernode_weight[hn_mapping[hn]] = nodeWeight(hn);
        }
      });
    }, [&] {
      if ( already_cut ) {
        extracted_block.already_cut.resize(num_hyperedges);
        const vec<uint8_t>& already_cut_hes = *already_cut;
        doParallelForAllEdges([&](const HyperedgeID he) {
          if ( he_mapping[he] != kInvalidHyperedge ) {
            ASSERT(he_mapping[he] < num_hyperedges);
            extracted_block.already_cut[he_mapping[he]] = already_cut_hes[he];
          }
        });
      }
    });

    // Construct hypergraph
    extracted_block.hg = HypergraphFactory::construct(num_hypernodes, num_hyperedges,
      edge_vector, hyperedge_weight.data(), hypernode_weight.data(), stable_construction_of_incident_edges);

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& hn) {
      if ( partID(hn) == block ) {
        const HypernodeID extracted_hn = hn_mapping[hn];
        extracted_block.hg.setCommunityID(extracted_hn, _hg->communityID(hn));
      }
    });
    return extracted_block;
  }

  // ! Extracts all blocks of the partition (from block 0 to block k).
  // ! This function has running time linear in the size of the original hypergraph
  // ! and should be used instead of extract(...) when more than two blocks should be extracted.
  std::pair<vec<ExtractedBlock>, vec<HypernodeID>> extractAllBlocks(const PartitionID k,
                                                                    const vec<uint8_t>* already_cut,
                                                                    const bool cut_net_splitting,
                                                                    const bool stable_construction_of_incident_edges) {
    /// [debug] std::cerr << "PartitionedHypergraph::extractAllBlocks(k, already_cut, cut_net_splitting, stable_construction_of_incident_edges)" << std::endl;
    ASSERT(k <= _k);

    vec<HypernodeID> hn_mapping(_hg->initialNumNodes(), kInvalidHypernode);
    vec<parallel::AtomicWrapper<HypernodeID>> nodes_cnt(
      k, parallel::AtomicWrapper<HypernodeID>(0));
    vec<vec<HyperedgeID>> hes2block(k);

    if ( stable_construction_of_incident_edges ) {
      // Stable construction for deterministic behavior requires
      // to determine node and edge IDs sequentially
      tbb::parallel_invoke([&] {
        // Compactify node IDs
        for ( const HypernodeID& hn : nodes() ) {
          const PartitionID block = partID(hn);
          if ( block < k ) {
            hn_mapping[hn] = nodes_cnt[block]++;
          }
        }
      }, [&] {
        // Get hyperedges contained in each block
        for ( const HyperedgeID& he : edges() ) {
          for ( const PartitionID& block : connectivitySet(he) ) {
            if ( pinCountInPart(he, block) > 1 &&
                (cut_net_splitting || connectivity(he) == 1) ) {
              hes2block[block].push_back(he);
            }
          }
        }
      });
    } else {
      vec<ds::StreamingVector<HyperedgeID>> hes2block_stream(k);
      tbb::parallel_invoke([&] {
        // Compactify node IDs
        doParallelForAllNodes([&](const HypernodeID& hn) {
          const PartitionID block = partID(hn);
          if ( block < k ) {
            hn_mapping[hn] = nodes_cnt[block]++;
          }
        });
      }, [&] {
        // Get hyperedges contained in each block
        doParallelForAllEdges([&](const HyperedgeID& he) {
          for ( const PartitionID& block : connectivitySet(he) ) {
            if ( pinCountInPart(he, block) > 1 &&
                (cut_net_splitting || connectivity(he) == 1) ) {
              hes2block_stream[block].stream(he);
            }
          }
        });
        // Copy hyperedges of a block into one vector
        tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
          hes2block[p] = hes2block_stream[p].copy_parallel();
        });
      });
    }

    // Extract plain hypergraph data for corresponding block
    using HyperedgeVector = vec<vec<HypernodeID>>;
    vec<ExtractedBlock> extracted_blocks(k);
    vec<HyperedgeVector> edge_vector(k);
    vec<vec<HyperedgeWeight>> he_weight(k);
    vec<vec<HypernodeWeight>> hn_weight(k);
    // Allocate auxilliary graph data structures
    tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
      const HypernodeID num_nodes = nodes_cnt[p];
      const HyperedgeID num_edges = hes2block[p].size();
      tbb::parallel_invoke([&] {
        edge_vector[p].resize(num_edges);
      }, [&] {
        he_weight[p].resize(num_edges);
      }, [&] {
        hn_weight[p].resize(num_nodes);
      }, [&] {
        if ( already_cut ) {
          extracted_blocks[p].already_cut.resize(num_edges);
        }
      });
    });

    // Write blocks to auxilliary graph data structure
    tbb::parallel_invoke([&] {
      tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
        tbb::parallel_for(UL(0), hes2block[p].size(), [&, p](const size_t i) {
          const HyperedgeID he = hes2block[p][i];
          he_weight[p][i] = edgeWeight(he);
          for ( const HypernodeID& pin : pins(he) ) {
            if ( partID(pin) == p ) {
              edge_vector[p][i].push_back(hn_mapping[pin]);
            }
          }
        });
      });
    }, [&] {
      doParallelForAllNodes([&](const HypernodeID& hn) {
        const PartitionID block = partID(hn);
        const HypernodeID mapped_hn = hn_mapping[hn];
        if ( block < k ) {
          ASSERT(UL(mapped_hn) < hn_weight[block].size());
          hn_weight[block][mapped_hn] = nodeWeight(hn);
        }
      });
    }, [&] {
      if ( already_cut ) {
        const vec<uint8_t>& already_cut_hes = *already_cut;
        tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
          tbb::parallel_for(UL(0), hes2block[p].size(), [&, p](const size_t i) {
            extracted_blocks[p].already_cut[i] = already_cut_hes[hes2block[p][i]];
          });
        });
      }
    });

    tbb::parallel_for(static_cast<PartitionID>(0), k, [&](const PartitionID p) {
      const HypernodeID num_nodes = nodes_cnt[p];
      const HyperedgeID num_hyperedges = hes2block[p].size();
      extracted_blocks[p].hg = HypergraphFactory::construct(num_nodes, num_hyperedges,
        edge_vector[p], he_weight[p].data(), hn_weight[p].data(),
        stable_construction_of_incident_edges);
    });

    // Set community ids
    doParallelForAllNodes([&](const HypernodeID& hn) {
      const PartitionID block = partID(hn);
      if ( block < k ) {
        extracted_blocks[block].hg.setCommunityID(hn_mapping[hn], _hg->communityID(hn));
      }
    });

    parallel::parallel_free(edge_vector);
    parallel::parallel_free(hn_weight, he_weight);

    return std::make_pair(std::move(extracted_blocks), std::move(hn_mapping));
  }

  void freeInternalData() {
    /// [debug] std::cerr << "PartitionedHypergraph::freeInternalData()" << std::endl;
    if ( _k > 0 ) {
      tbb::parallel_invoke( [&] {
        parallel::parallel_free(_part_ids, _pin_count_update_ownership);
      }, [&] {
        _con_info.freeInternalData();
      } );
    }
    _k = 0;
  }

 private:
  void applyPartWeightUpdates(vec<HypernodeWeight>& part_weight_deltas) {
    /// [debug] std::cerr << "PartitionedHypergraph::applyPartWeightUpdates(part_weight_deltas)" << std::endl;
    for (PartitionID p = 0; p < _k; ++p) {
      _part_weights[p].fetch_add(part_weight_deltas[p], std::memory_order_relaxed);
    }
  }

  void applyPartVolumeUpdates(vec<HyperedgeWeight>& part_volume_deltas) {
    /// [debug] std::cerr << "PartitionedHypergraph::applyPartVolumeUpdates(part_volume_deltas)" << std::endl;
    for (PartitionID p = 0; p < _k; ++p) {
      _part_volumes[p].fetch_add(part_volume_deltas[p], std::memory_order_relaxed);
    }
  }

  void applyPartCutWeightUpdates(vec<HyperedgeWeight>& part_cut_weight_deltas) {
    /// [debug] std::cerr << "PartitionedHypergraph::applyPartCutWeightUpdates(part_cut_weight_deltas)" << std::endl;
    for (PartitionID p = 0; p < _k; ++p) {
      _part_cut_weights[p].fetch_add(part_cut_weight_deltas[p], std::memory_order_relaxed);
    }
  }

  void initializeBlockWeights() {
    /// [debug] std::cerr << "PartitionedHypergraph::initializeBlockWeights()" << std::endl;
    auto accumulate = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HypernodeWeight> pws(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( nodeIsEnabled(u) ) {
          const PartitionID pu = partID( u );
          const HypernodeWeight wu = nodeWeight( u );
          pws[pu] += wu;
        }
      }
      applyPartWeightUpdates(pws);
    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  void initializeBlockVolumes() {
    /// [debug] std::cerr << "PartitionedHypergraph::initializeBlockVolumes()" << std::endl;
    auto accumulate = [&](tbb::blocked_range<HypernodeID>& r) {
      vec<HyperedgeWeight> pvs(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HypernodeID u = r.begin(); u < r.end(); ++u) {
        if ( nodeIsEnabled(u) ) {
          const PartitionID pu = partID( u );
          const HyperedgeWeight uWDeg = nodeWeightedDegree( u ); 
          // TODO: make sure, disabled edges are not a problem for weighted degree
          pvs[pu] += uWDeg;
        }
      }
      applyPartVolumeUpdates(pvs);
    };

    tbb::parallel_for(tbb::blocked_range<HypernodeID>(HypernodeID(0), initialNumNodes()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  void initializePinCountInPart() {
    /// [debug] std::cerr << "PartitionedHypergraph::initializePinCountInPart()" << std::endl;
    tls_enumerable_thread_specific< vec<HypernodeID> > ets_pin_count_in_part(_k, 0);

    auto assign = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HypernodeID>& pin_counts = ets_pin_count_in_part.local();
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {
        if ( edgeIsEnabled(he) ) {
          for (const HypernodeID& pin : pins(he)) {
            ++pin_counts[partID(pin)];
          }

          for (PartitionID p = 0; p < _k; ++p) {
            ASSERT(pinCountInPart(he, p) == 0);
            if (pin_counts[p] > 0) {
              _con_info.addBlock(he, p);
              _con_info.setPinCountInPart(he, p, pin_counts[p]);
            }
            pin_counts[p] = 0;
          }
        }
      }
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HyperedgeID(0), initialNumEdges()), assign);
  }

  // ! Uses connectivity, connectivitySet
  // ! => use only after initializePinCountInPart !!!
  void initializeBlockCutWeights() {
    /// [debug] std::cerr << "PartitionedHypergraph::initializeBlockCutWeights()" << std::endl;
    auto accumulate = [&](tbb::blocked_range<HyperedgeID>& r) {
      vec<HyperedgeWeight> pcws(_k, 0);  // this is not enumerable_thread_specific because of the static partitioner
      for (HyperedgeID he = r.begin(); he < r.end(); ++he) {
        if ( edgeIsEnabled(he) && connectivity(he) > 1 ) {
          for ( PartitionID partId : connectivitySet(he) ) {
            pcws[partId] += edgeWeight(he);
          }
        }
      }
      applyPartCutWeightUpdates(pcws);
    };

    tbb::parallel_for(tbb::blocked_range<HyperedgeID>(HyperedgeID(0), initialNumEdges()),
                      accumulate,
                      tbb::static_partitioner()
    );
  }

  HypernodeID pinCountInPartRecomputed(const HyperedgeID e, PartitionID p) const {
    HypernodeID pcip = 0;
    for (HypernodeID u : pins(e)) {
      if (partID(u) == p) {
        pcip++;
      }
    }
    return pcip;
  }

  // ! Updates pin count in part using a spinlock.
  // ! Also updates _part_cut_weights (after my adjustments)
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void updatePinCountOfHyperedge(const HyperedgeID he,
                                                                    const PartitionID from,
                                                                    const PartitionID to,
                                                                    SynchronizedEdgeUpdate& sync_update,
                                                                    const DeltaFunction& delta_func,
                                                                    const NotificationFunc& notify_func) {
    /// [debug] std::cerr << "PartitionedHypergraph::updatePinCountOfHyperedge(he, from, to, sync_update, delta_func, notify_func)" << std::endl;
    ASSERT(he < _pin_count_update_ownership.size());
    // TODO: change sync_update to contain info for calculating concuctance gain (1st and 2nd max cond cuts; total val; ...)
    //       Problem: conductance gain depends on both whole parts, not just the edge
    sync_update.he = he;
    sync_update.edge_weight = edgeWeight(he);
    sync_update.edge_size = edgeSize(he);
    _pin_count_update_ownership[he].lock();
    notify_func(sync_update);
    const HypernodeID old_pins_in_from_part = pinCountInPart(he, from);
    sync_update.pin_count_in_from_part_after = decrementPinCountOfBlock(he, from);
    sync_update.pin_count_in_to_part_after = incrementPinCountOfBlock(he, to);
    sync_update.connectivity_set_after = hasTargetGraph() ? &deepCopyOfConnectivitySet(he) : nullptr;
    sync_update.pin_counts_after = hasTargetGraph() ? &_con_info.pinCountSnapshot(he) : nullptr;
    const HypernodeID new_pins_in_to_part = pinCountInPart(he, to);
    // update _part_cut_weights for "from" part
    if (HypernodeID(1) == old_pins_in_from_part && old_pins_in_from_part < edgeSize(he)) {
      // he was a cutting edge for part "from", but not anymore
      _part_cut_weights[from].fetch_sub(edgeWeight(he), std::memory_order_relaxed);
    } else if (HypernodeID(1) < old_pins_in_from_part && old_pins_in_from_part == edgeSize(he)) {
      // he was not a cutting edge for part "from", but now is
      _part_cut_weights[from].fetch_add(edgeWeight(he), std::memory_order_relaxed);
    }
    // update _part_cut_weights for "to" part
    if (HypernodeID(1) == new_pins_in_to_part && new_pins_in_to_part < edgeSize(he)) {
      // he was not a cutting edge for part to, but now is
      _part_cut_weights[to].fetch_add(edgeWeight(he), std::memory_order_relaxed);
    } else if (HypernodeID(1) < new_pins_in_to_part && new_pins_in_to_part == edgeSize(he)) {
      // he was a cutting edge for part "to", but not anymore
      _part_cut_weights[to].fetch_sub(edgeWeight(he), std::memory_order_relaxed);
    }
    // for all other parts, _part_cut_weights remains the same
    _pin_count_update_ownership[he].unlock();
    delta_func(sync_update);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID decrementPinCountOfBlock(const HyperedgeID e, const PartitionID p) {
    /// [debug] std::cerr << "PartitionedHypergraph::decrementPinCountOfBlock(e, p)" << std::endl;
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _con_info.decrementPinCountInPart(e, p);
    if ( pin_count_after == 0 ) {
      _con_info.removeBlock(e, p);
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HypernodeID incrementPinCountOfBlock(const HyperedgeID e, const PartitionID p) {
    /// [debug] std::cerr << "PartitionedHypergraph::incrementPinCountOfBlock(e, p)" << std::endl;
    ASSERT(e < _hg->initialNumEdges(), "Hyperedge" << e << "does not exist");
    ASSERT(edgeIsEnabled(e), "Hyperedge" << e << "is disabled");
    ASSERT(p != kInvalidPartition && p < _k);
    const HypernodeID pin_count_after = _con_info.incrementPinCountInPart(e, p);
    if ( pin_count_after == 1 ) {
      _con_info.addBlock(e, p);
    }
    return pin_count_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight decrementCutWeightOfBlock(const PartitionID p, const HyperedgeWeight w) {
    /// [debug] std::cerr << "PartitionedHypergraph::decrementCutWeightOfBlock(p, w)" << std::endl;
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeWeight cut_weight_after = _part_cut_weights[p].fetch_sub(w, std::memory_order_relaxed);
    return cut_weight_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight incrementCutWeightOfBlock(const PartitionID p, const HyperedgeWeight w) {
    /// [debug] std::cerr << "PartitionedHypergraph::incrementCutWeightOfBlock(p, w)" << std::endl;
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeWeight cut_weight_after = _part_cut_weights[p].fetch_add(w, std::memory_order_relaxed);
    return cut_weight_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight decrementVolumeOfBlock(const PartitionID p, const HyperedgeWeight w) {
    /// [debug] std::cerr << "PartitionedHypergraph::decrementVolumeOfBlock(p, w)" << std::endl;
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeWeight volume_after = _part_volumes[p].fetch_sub(w, std::memory_order_relaxed);
    return volume_after;
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE
  HyperedgeWeight incrementVolumeOfBlock(const PartitionID p, const HyperedgeWeight w) {
    /// [debug] std::cerr << "PartitionedHypergraph::incrementVolumeOfBlock(p, w)" << std::endl;
    ASSERT(p != kInvalidPartition && p < _k);
    const HyperedgeWeight volume_after = _part_volumes[p].fetch_add(w, std::memory_order_relaxed);
    return volume_after;
  }

  // ! Number of nodes of the top level hypergraph
  HypernodeID _input_num_nodes = 0;

  // ! Number of hyperedges of the top level hypergraph
  HyperedgeID _input_num_edges = 0;

  // ! Number of blocks
  PartitionID _k = 0;

  // ! Underlying hypergraph
  Hypergraph* _hg = nullptr;

  // ! Target graph on which this hypergraph is mapped
  const TargetGraph* _target_graph;

  // ! Conductance PQ (needs to be enabled)
  ConductancePriorityQueue<Self> _conductance_pq;

  // ! Weight and information for all blocks.
  vec< CAtomic<HypernodeWeight> > _part_weights;
  
  // ! Volume and information for all blocks.
  vec< CAtomic<HypernodeWeight> > _part_volumes;
  
  // ! Sum of weights of cutting edges for all blocks.
  vec< CAtomic<HypernodeWeight> > _part_cut_weights;

  // ! Current block IDs of the vertices
  Array< PartitionID > _part_ids;

  // ! Stores the pin count values and connectivity sets
  ConnectivityInformation _con_info;

  // ! In order to update the pin count of a hyperedge thread-safe, a thread must acquire
  // ! the ownership of a hyperedge via a CAS operation.
  Array<SpinLock> _pin_count_update_ownership;
};

} // namespace ds
} // namespace mt_kahypar
