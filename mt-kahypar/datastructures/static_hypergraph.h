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


#include <tbb/parallel_for.h>

#include "include/mtkahypartypes.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/fixed_vertex_support.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/partition/context_enum_classes.h"
#include "mt-kahypar/utils/memory_tree.h"
#include "mt-kahypar/utils/range.h"
#include "mt-kahypar/utils/exception.h"

namespace mt_kahypar {
namespace ds {

// Forward
class StaticHypergraphFactory;
template <typename Hypergraph,
          typename ConnectivityInformation>
class PartitionedHypergraph;

class StaticHypergraph {

  static constexpr bool enable_heavy_assert = false;

  // During contractions we temporary memcpy all incident nets of a collapsed
  // vertex to consecutive range in a temporary incident nets structure.
  // Afterwards, we sort that range and remove duplicates. However, it turned
  // out that this become a major sequential bottleneck in presence of high
  // degree vertices. Therefore, all vertices with temporary degree greater
  // than this threshold are contracted with a special procedure.
  static constexpr HyperedgeID HIGH_DEGREE_CONTRACTION_THRESHOLD = ID(500000);

  static_assert(std::is_unsigned<HypernodeID>::value, "Hypernode ID must be unsigned");
  static_assert(std::is_unsigned<HyperedgeID>::value, "Hyperedge ID must be unsigned");

  using AtomicHypernodeID = parallel::IntegralAtomicWrapper<HypernodeID>;
  using AtomicHypernodeWeight = parallel::IntegralAtomicWrapper<HypernodeWeight>;
  using UncontractionFunction = std::function<void (const HypernodeID, const HypernodeID, const HyperedgeID)>;
  #define NOOP_BATCH_FUNC [] (const HypernodeID, const HypernodeID, const HyperedgeID) { }

  /**
   * Represents a hypernode of the hypergraph and contains all information
   * associated with a vertex.
   */
  class Hypernode {
   public:
    using IDType = HypernodeID;

    Hypernode() :
      _begin(0),
      _size(0),
      _weight(1),
      _valid(false) { }

    Hypernode(const bool valid) :
      _begin(0),
      _size(0),
      _weight(1),
      _valid(valid) { }

    // Sentinel Constructor
    Hypernode(const size_t begin) :
      _begin(begin),
      _size(0),
      _weight(1),
      _valid(false) { }

    bool isDisabled() const {
      return _valid == false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    // ! Returns the index of the first element in _incident_nets
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incident_nets to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incident_nets
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    HyperedgeWeight weight() const {
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

   private:
    // ! Index of the first element in _incident_nets
    size_t _begin;
    // ! Number of incident nets
    size_t _size;
    // ! Hypernode weight
    HypernodeWeight _weight;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /**
   * Represents a hyperedge of the hypergraph and contains all information
   * associated with a net (except connectivity information).
   */
  class Hyperedge {
   public:
    using IDType = HyperedgeID;

    Hyperedge() :
      _begin(0),
      _size(0),
      _original_size(0),
      _weight(1),
      _valid(false) { }

    // Sentinel Constructor
    Hyperedge(const size_t begin) :
      _begin(begin),
      _size(0),
      _original_size(0),
      _weight(1),
      _valid(false) { }

    // ! Disables the hypernode/hyperedge. Disable hypernodes/hyperedges will be skipped
    // ! when iterating over the set of all nodes/edges.
    void disable() {
      ASSERT(!isDisabled());
      _valid = false;
    }

    void enable() {
      ASSERT(isDisabled());
      _valid = true;
    }

    bool isDisabled() const {
      return _valid == false;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstEntry() const {
      return _begin;
    }

    // ! Sets the index of the first element in _incidence_array to begin
    void setFirstEntry(size_t begin) {
      ASSERT(!isDisabled());
      _begin = begin;
    }

    // ! Returns the index of the first element in _incidence_array
    size_t firstInvalidEntry() const {
      return _begin + _size;
    }

    size_t size() const {
      ASSERT(!isDisabled());
      return _size;
    }

    void setSize(size_t size) {
      ASSERT(!isDisabled());
      _size = size;
    }

    size_t originalSize() const {
      ASSERT(!isDisabled());
      return _original_size;
    }

    void setOriginalSize(size_t original_size) {
      ASSERT(!isDisabled());
      _original_size = original_size;
    }

    HyperedgeWeight weight() const {
      ASSERT(!isDisabled());
      return _weight;
    }

    void setWeight(HyperedgeWeight weight) {
      ASSERT(!isDisabled());
      _weight = weight;
    }

    bool operator== (const Hyperedge& rhs) const {
      return _begin == rhs._begin && _size == rhs._size && _weight == rhs._weight;
    }

    bool operator!= (const Hyperedge& rhs) const {
      return _begin != rhs._begin || _size != rhs._size || _weight != rhs._weight;
    }

   private:
    // ! Index of the first element in _incidence_array
    size_t _begin;
    // ! Number of pins
    size_t _size;
    // ! Number of pins **at the moment of snapshot**
    size_t _original_size;
    // ! hyperedge weight
    HyperedgeWeight _weight;
    // ! Flag indicating whether or not the element is active.
    bool _valid;
  };

  /*!
   * Iterator for HypergraphElements (Hypernodes/Hyperedges)
   *
   * The iterator is used in for-each loops over all hypernodes/hyperedges.
   * In order to support iteration over coarsened hypergraphs, this iterator
   * skips over HypergraphElements marked as invalid.
   * Iterating over the set of vertices \f$V\f$ therefore is linear in the
   * size \f$|V|\f$ of the original hypergraph - even if it has been coarsened
   * to much smaller size. The same also holds true for for-each loops over
   * the set of hyperedges.
   *
   * In order to be as generic as possible, the iterator does not expose the
   * internal Hypernode/Hyperedge representations. Instead only handles to
   * the respective elements are returned, i.e. the IDs of the corresponding
   * hypernodes/hyperedges.
   *
   */
  template <typename ElementType>
  class HypergraphElementIterator {
   public:
    using IDType = typename ElementType::IDType;
    using iterator_category = std::forward_iterator_tag;
    using value_type = IDType;
    using reference = IDType&;
    using pointer = const IDType*;
    using difference_type = std::ptrdiff_t;

    /*!
     * Construct a HypergraphElementIterator
     * See GenericHypergraph::nodes() or GenericHypergraph::edges() for usage.
     *
     * If start_element is invalid, the iterator advances to the first valid
     * element.
     *
     * \param start_element A pointer to the starting position
     * \param id The index of the element the pointer points to
     * \param max_id The maximum index allowed
     */
    HypergraphElementIterator(const ElementType* start_element, IDType id, IDType max_id) :
      _id(id),
      _max_id(max_id),
      _element(start_element) {
      if (_id != _max_id && _element->isDisabled()) {
        operator++ ();
      }
    }

    // ! Returns the id of the element the iterator currently points to.
    IDType operator* () const {
      return _id;
    }

    // ! Prefix increment. The iterator advances to the next valid element.
    HypergraphElementIterator & operator++ () {
      ASSERT(_id < _max_id);
      do {
        ++_id;
        ++_element;
      } while (_id < _max_id && _element->isDisabled());
      return *this;
    }

    // ! Postfix increment. The iterator advances to the next valid element.
    HypergraphElementIterator operator++ (int) {
      HypergraphElementIterator copy = *this;
      operator++ ();
      return copy;
    }

    bool operator!= (const HypergraphElementIterator& rhs) {
      return _id != rhs._id;
    }

    bool operator== (const HypergraphElementIterator& rhs) {
      return _id == rhs._id;
    }

   private:
    // Handle to the HypergraphElement the iterator currently points to
    IDType _id = 0;
    // Maximum allowed index
    IDType _max_id = 0;
    // HypergraphElement the iterator currently points to
    const ElementType* _element = nullptr;
  };

  static_assert(std::is_trivially_copyable<Hypernode>::value, "Hypernode is not trivially copyable");
  static_assert(std::is_trivially_copyable<Hyperedge>::value, "Hyperedge is not trivially copyable");

  using IncidenceArray = Array<HypernodeID>;
  using IncidentNets = Array<HyperedgeID>;

  // ! Contains buffers that are needed during multilevel contractions.
  // ! Struct is allocated on top level hypergraph and passed to each contracted
  // ! hypergraph such that memory can be reused in consecutive contractions.
  struct TmpContractionBuffer {
    explicit TmpContractionBuffer(const HypernodeID num_hypernodes,
                                  const HyperedgeID num_hyperedges,
                                  const HyperedgeID num_pins) {
      tbb::parallel_invoke([&] {
        mapping.resize("Coarsening", "mapping", num_hypernodes);
      }, [&] {
        tmp_hypernodes.resize("Coarsening", "tmp_hypernodes", num_hypernodes);
      }, [&] {
        tmp_incident_nets.resize("Coarsening", "tmp_incident_nets", num_pins);
      }, [&] {
        tmp_num_incident_nets.resize("Coarsening", "tmp_num_incident_nets", num_hypernodes);
      }, [&] {
        hn_weights.resize("Coarsening", "hn_weights", num_hypernodes);
      }, [&] {
        tmp_hyperedges.resize("Coarsening", "tmp_hyperedges", num_hyperedges);
      }, [&] {
        tmp_incidence_array.resize("Coarsening", "tmp_incidence_array", num_pins);
      }, [&] {
        he_sizes.resize("Coarsening", "he_sizes", num_hyperedges);
      }, [&] {
        valid_hyperedges.resize("Coarsening", "valid_hyperedges", num_hyperedges);
      });
    }

    Array<size_t> mapping;
    Array<Hypernode> tmp_hypernodes;
    IncidentNets tmp_incident_nets;
    Array<parallel::IntegralAtomicWrapper<size_t>> tmp_num_incident_nets;
    Array<parallel::IntegralAtomicWrapper<HypernodeWeight>> hn_weights;
    Array<Hyperedge> tmp_hyperedges;
    IncidenceArray tmp_incidence_array;
    Array<size_t> he_sizes;
    Array<size_t> valid_hyperedges;
  };

 public:
  static constexpr bool is_graph = false;
  static constexpr bool is_static_hypergraph = true;
  static constexpr bool is_partitioned = false;
  static constexpr size_t SIZE_OF_HYPERNODE = sizeof(Hypernode);
  static constexpr size_t SIZE_OF_HYPEREDGE = sizeof(Hyperedge);
  static constexpr mt_kahypar_hypergraph_type_t TYPE = STATIC_HYPERGRAPH;

  // ! Factory
  using Factory = StaticHypergraphFactory;
  // ! Iterator to iterate over the hypernodes
  using HypernodeIterator = HypergraphElementIterator<const Hypernode>;
  // ! Iterator to iterate over the hyperedges
  using HyperedgeIterator = HypergraphElementIterator<const Hyperedge>;
  // ! Iterator to iterate over the pins of a hyperedge
  using IncidenceIterator = typename IncidenceArray::const_iterator;
  // ! Iterator to iterate over the incident nets of a hypernode
  using IncidentNetsIterator = typename IncidentNets::const_iterator;

  struct ParallelHyperedge {
    HyperedgeID removed_hyperedge;
    HyperedgeID representative;
  };

  explicit StaticHypergraph() :
    _num_hypernodes(0),
    _num_removed_hypernodes(0),
    _removed_degree_zero_hn_weight(0),
    _num_hyperedges(0),
    _num_removed_hyperedges(0),
    _max_edge_size(0),
    _num_pins(0),
    _total_degree(0),
    _total_weight(0),
    _total_volume(0),
    _original_total_volume(0),
    _hypernodes(),
    _incident_nets(),
    _hyperedges(),
    _incidence_array(),
    _weighted_degrees(),
    _original_weighted_degrees(),
    _community_ids(0),
    _fixed_vertices(),
    _tmp_contraction_buffer(nullptr),
    _beta(1, 0.0), _gamma(1, 0.0), _omega(1, {0.0, 0.0})
    { }

  StaticHypergraph(const StaticHypergraph&) = delete;
  StaticHypergraph & operator= (const StaticHypergraph &) = delete;

  StaticHypergraph(StaticHypergraph&& other) :
    _num_hypernodes(other._num_hypernodes),
    _num_removed_hypernodes(other._num_removed_hypernodes),
    _removed_degree_zero_hn_weight(other._removed_degree_zero_hn_weight),
    _num_hyperedges(other._num_hyperedges),
    _num_removed_hyperedges(other._num_removed_hyperedges),
    _max_edge_size(other._max_edge_size),
    _num_pins(other._num_pins),
    _total_degree(other._total_degree),
    _total_weight(other._total_weight),
    _total_volume(other._total_volume),
    _original_total_volume(other._original_total_volume),
    _hypernodes(std::move(other._hypernodes)),
    _incident_nets(std::move(other._incident_nets)),
    _hyperedges(std::move(other._hyperedges)),
    _incidence_array(std::move(other._incidence_array)),
    _weighted_degrees(std::move(other._weighted_degrees)),
    _original_weighted_degrees(std::move(other._original_weighted_degrees)),
    _community_ids(std::move(other._community_ids)),
    _fixed_vertices(std::move(other._fixed_vertices)),
    _tmp_contraction_buffer(std::move(other._tmp_contraction_buffer)),
    _enable_collective_sync_update(other._enable_collective_sync_update),
    _disable_single_pin_nets_removal(other._disable_single_pin_nets_removal),
    _beta(std::move(other._beta)), _gamma(std::move(other._gamma)),
        _omega(std::move(other._omega)) 
  {
    _fixed_vertices.setHypergraph(this);
    other._tmp_contraction_buffer = nullptr;
  }

  StaticHypergraph & operator= (StaticHypergraph&& other) {
    _num_hypernodes = other._num_hypernodes;
    _num_removed_hypernodes = other._num_removed_hypernodes;
    _removed_degree_zero_hn_weight = other._removed_degree_zero_hn_weight;
    _num_hyperedges = other._num_hyperedges;
    _num_removed_hyperedges = other._num_removed_hyperedges;
    _max_edge_size = other._max_edge_size;
    _num_pins = other._num_pins;
    _total_degree = other._total_degree;
    _total_weight = other._total_weight;
    _total_volume = other._total_volume;
    _original_total_volume = other._original_total_volume;
    _hypernodes = std::move(other._hypernodes);
    _incident_nets = std::move(other._incident_nets);
    _hyperedges = std::move(other._hyperedges);
    _incidence_array = std::move(other._incidence_array);
    _weighted_degrees = std::move(other._weighted_degrees);
    _original_weighted_degrees = std::move(other._original_weighted_degrees);
    _community_ids = std::move(other._community_ids);
    _fixed_vertices = std::move(other._fixed_vertices);
    _fixed_vertices.setHypergraph(this);
    _tmp_contraction_buffer = std::move(other._tmp_contraction_buffer);
    _disable_single_pin_nets_removal = other._disable_single_pin_nets_removal;
    _enable_collective_sync_update = other._enable_collective_sync_update;
    other._tmp_contraction_buffer = nullptr;
    _beta = std::move(other._beta);
    _gamma = std::move(other._gamma);
    _omega = std::move(other._omega);
    return *this;
  }

  ~StaticHypergraph() {
    if ( _tmp_contraction_buffer ) {
      delete(_tmp_contraction_buffer);
      _tmp_contraction_buffer = nullptr;
    }
    freeInternalData();
  }

  // ####################### General Hypergraph Stats #######################

  // ! Initial number of hypernodes
  HypernodeID initialNumNodes() const {
    return _num_hypernodes;
  }

  // ! Number of removed hypernodes
  HypernodeID numRemovedHypernodes() const {
    return _num_removed_hypernodes;
  }

  // ! Weight of removed degree zero vertics
  HypernodeWeight weightOfRemovedDegreeZeroVertices() const {
    return _removed_degree_zero_hn_weight;
  }

  // ! Initial number of hyperedges
  HyperedgeID initialNumEdges() const {
    return _num_hyperedges;
  }

  // ! Number of removed hyperedges
  HyperedgeID numRemovedHyperedges() const {
    return _num_removed_hyperedges;
  }

  // ! Set the number of removed hyperedges
  void setNumRemovedHyperedges(const HyperedgeID num_removed_hyperedges) {
    _num_removed_hyperedges = num_removed_hyperedges;
  }

  // ! Initial number of pins
  HypernodeID initialNumPins() const {
    return _num_pins;
  }

  // ! Initial sum of the degree of all vertices
  HypernodeID initialTotalVertexDegree() const {
    return _total_degree;
  }

  // ! Total weight of hypergraph
  HypernodeWeight totalWeight() const {
    return _total_weight;
  }

  // ! Total volume of hypergraph
  HypergraphVolume totalVolume() const {
    return _total_volume;
  }

  // ! Original total volume of hypergraph
  HypergraphVolume originalTotalVolume() const {
    return _original_total_volume;
  }

  // ! Computes the total node weight of the hypergraph
  void computeAndSetTotalNodeWeight(parallel_tag_t);

  // ! Computes the total volume of the hypergraph
  void computeAndSetTotalVolume(parallel_tag_t);

  // ####################### Iterators #######################

  // ! Iterates in parallel over all active nodes and calls function f
  // ! for each vertex
  template<typename F>
  void doParallelForAllNodes(const F& f) const {
    tbb::parallel_for(ID(0), _num_hypernodes, [&](const HypernodeID& hn) {
      if ( nodeIsEnabled(hn) ) {
        f(hn);
      }
    });
  }

  // ! Iterates in parallel over all active edges and calls function f
  // ! for each net
  template<typename F>
  void doParallelForAllEdges(const F& f) const {
    tbb::parallel_for(ID(0), _num_hyperedges, [&](const HyperedgeID& he) {
      if ( edgeIsEnabled(he) ) {
        f(he);
      }
    });
  }

  // ! Returns a range of the active nodes of the hypergraph
  IteratorRange<HypernodeIterator> nodes() const {
    return IteratorRange<HypernodeIterator>(
      HypernodeIterator(_hypernodes.data(), ID(0), _num_hypernodes),
      HypernodeIterator(_hypernodes.data() + _num_hypernodes, _num_hypernodes, _num_hypernodes));
  }

  // ! Returns a range of the active edges of the hypergraph
  IteratorRange<HyperedgeIterator> edges() const {
    return IteratorRange<HyperedgeIterator>(
      HyperedgeIterator(_hyperedges.data(), ID(0), _num_hyperedges),
      HyperedgeIterator(_hyperedges.data() + _num_hyperedges, _num_hyperedges, _num_hyperedges));
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetsIterator> incidentEdges(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    const Hypernode& hn = hypernode(u);
    return IteratorRange<IncidentNetsIterator>(
      _incident_nets.cbegin() + hn.firstEntry(),
      _incident_nets.cbegin() + hn.firstInvalidEntry());
  }

  // ! Returns a range to loop over the pins of hyperedge e.
  IteratorRange<IncidenceIterator> pins(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    const Hyperedge& he = hyperedge(e);
    return IteratorRange<IncidenceIterator>(
      _incidence_array.cbegin() + he.firstEntry(),
      _incidence_array.cbegin() + he.firstInvalidEntry());
  }

    // ####################### Hypernode Information #######################

  // ! Weight of a vertex
  HypernodeWeight nodeWeight(const HypernodeID u) const {
    return hypernode(u).weight();
  }

  // ! Sets the weight of a vertex
  void setNodeWeight(const HypernodeID u, const HypernodeWeight weight) {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).setWeight(weight);
  }

  // ! Degree of a hypernode
  HyperedgeID nodeDegree(const HypernodeID u) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    return hypernode(u).size();
  }

  // ! Weighted degree of a hypernode
  HypergraphVolume nodeWeightedDegree(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _weighted_degrees[u];
  }

  // ! Decrease weighted degree of a hypernode
  // ! (Not supported)
  void decreaseNodeWeightedDegree(const HypernodeID u, const HypergraphVolume w) const {
    unused(u);
    unused(w);
    throw UnsupportedOperationException(
      "decreaseNodeWeightedDegree(u, w) is not supported in static hypergraph");
  }

  // ! Original weighted degree of a hypernode
  // ! (during dontractions and single-pin nets removal)
  HypergraphVolume nodeOriginalWeightedDegree(const HypernodeID u) const {
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return _original_weighted_degrees[u];
  }
  
  // ! Returns, whether a hypernode is enabled or not
  bool nodeIsEnabled(const HypernodeID u) const {
    return !hypernode(u).isDisabled();
  }

  // ! Enables a hypernode (must be disabled before)
  void enableHypernode(const HypernodeID u) {
    hypernode(u).enable();
  }

  // ! Disables a hypernode (must be enabled before)
  void disableHypernode(const HypernodeID u) {
    hypernode(u).disable();
  }

  // ! Removes a hypernode (must be enabled before)
  void removeHypernode(const HypernodeID u) {
    hypernode(u).disable();
    ++_num_removed_hypernodes;
  }

  // ! Removes a degree zero hypernode
  void removeDegreeZeroHypernode(const HypernodeID u) {
    ASSERT(nodeDegree(u) == 0);
    _removed_degree_zero_hn_weight += nodeWeight(u);
    removeHypernode(u);
  }

  // ! Restores a degree zero hypernode
  void restoreDegreeZeroHypernode(const HypernodeID u) {
    hypernode(u).enable();
    ASSERT(nodeDegree(u) == 0);
    _removed_degree_zero_hn_weight -= nodeWeight(u);
  }

  // ####################### Hyperedge Information #######################

  // ! Weight of a hyperedge
  HyperedgeWeight edgeWeight(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).weight();
  }

  // ! Sets the weight of a hyperedge
  void setEdgeWeight(const HyperedgeID e, const HyperedgeWeight weight) {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).setWeight(weight);
  }

  // ! Number of pins of a hyperedge
  HypernodeID edgeSize(const HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).size();
  }

  // ! Maximum size of a hyperedge
  HypernodeID maxEdgeSize() const {
    return _max_edge_size;
  }

  // ! Returns, whether a hyperedge is enabled or not
  bool edgeIsEnabled(const HyperedgeID e) const {
    return !hyperedge(e).isDisabled();
  }

  // ! Enables a hyperedge (must be disabled before)
  void enableHyperedge(const HyperedgeID e) {
    hyperedge(e).enable();
  }

  // ! Disabled a hyperedge (must be enabled before)
  void disableHyperedge(const HyperedgeID e) {
    hyperedge(e).disable();
  }

  // ! Community id which hypernode u is assigned to
  PartitionID communityID(const HypernodeID u) const {
    return _community_ids[u];
  }

  // ! Assign a community to a hypernode
  void setCommunityID(const HypernodeID u, const PartitionID community_id) {
    _community_ids[u] = community_id;
  }

  // ########################### Snapshots ############################
  // ----------- Snapshot edge sizes ----------------

  // ! Save current edge sizes as original edge sizes
  void snapshotOriginalEdgeSizes() {
    for (HyperedgeID e : edges()) {
      hyperedge(e).setOriginalSize(hyperedge(e).size());
    }
  }

  // ! Get the edge size at the moment of the last snapshot
  HypernodeID originalEdgeSize(HyperedgeID e) const {
    ASSERT(!hyperedge(e).isDisabled(), "Hyperedge" << e << "is disabled");
    return hyperedge(e).originalSize();
  }

  // ! Get the maximum edge size at the moment of the last snapshot
  HypernodeID originalMaxEdgeSize() const {
    return _original_max_edge_size;
  }

  // ------- Snapshot volumes and weighted degrees -------
private:
  // ! Save the current weighted degrees as original
  // ! (private as weighted degrees should be consistent with the total volume)
  void snapshotOriginalWeightedDegrees() {
    _original_weighted_degrees = _weighted_degrees;
  }
  // ! Save the current total volume as original
  // ! (private as total volume should be consistent with the weighted degrees)
  void snapshotOriginalTotalVolume() {
    _original_total_volume = totalVolume();
  }
public:
  // ! Save the current weighted degrees and total volume as original stats
  // ! (together for their consistency)
  void snapshotOriginalWeightedDegreesAndTotalVolume() {
    snapshotOriginalWeightedDegrees();
    snapshotOriginalTotalVolume();
  }

  // ######################## AON-Hypermodularity #######################
  
  // ! AON HyperModularity Clustering
  // ! true once the three vectors were filled at the finest level
  bool hasAON() const { return !_beta.empty(); }

  // ! Get \beta for AON HyperModularity Clustering
  // ! Constant-time access by edge size d (d ≥ 0, d < _beta.size())
  double beta(std::size_t d) const noexcept {
    ASSERT(0 <= d && d < _beta.size(),
           "d = " << d << " is out of bounds for beta vector of size " << _beta.size());
    return _beta[d];
    // return d < _beta .size() ? _beta [d]      : 0.0;
  }

  // ! Get \gamma for AON HyperModularity Clustering
  // ! Constant-time access by edge size d (d ≥ 0, d < _gamma.size())
  double gamma(std::size_t d) const noexcept {
    ASSERT(0 <= d && d < _gamma.size(),
           "d = " << d << " is out of bounds for gamma vector of size " << _gamma.size());
    return _gamma[d];
    // return d < _gamma.size() ? _gamma[d]      : 0.0;
  }

  // ! Get \omega_{d0} for AON HyperModularity Clustering
  // ! Constant-time access by edge size d (d ≥ 0, d < _omega.size())
  double omegaIn(std::size_t d) const noexcept {
    ASSERT(0 <= d && d < _omega.size(),
           "d = " << d << " is out of bounds for omega vector of size " << _omega.size());
    return _omega[d][0];
    // return d < _omega.size() ? _omega[d][0]   : 0.0;
  }

  // ! Get \omega_{d1} for AON HyperModularity Clustering
  // ! Constant-time access by edge size d (d ≥ 0, d < _omega.size())
  double omegaOut(std::size_t d) const noexcept {
    ASSERT(0 <= d && d < _omega.size(),
           "d = " << d << " is out of bounds for omega vector of size " << _omega.size());
    return _omega[d][1];
    // return d < _omega.size() ? _omega[d][1]   : 0.0;
  }

  // ──────────────────────────────────────────────────────────
  /// (Re)compute β, γ, ω for the **current** community assignment
  /// `_community_ids`.
  ///
  /// *   β[k]   = log(ω_in/ω_out)   for edge-size k (k ≥ 2)
  /// *   γ[k]   = ω_in − ω_out
  /// *   ω[k] = { ω_in , ω_out }
  /// 
  /// [Mariia: _beta = - \beta, _gamma = - \beta * \gamma 
  ///           with \beta, \gamma from the Hypermodularity article]
  ///
  /// After the call the three member vectors `_beta`, `_gamma`, `_omega` are
  /// filled and can be queried with beta(k), gamma(k), omegaIn/Out(k).
  inline void computeAONParameters(double eps = 1e-12) {
    const std::size_t dmax = static_cast<std::size_t>(_max_edge_size);
    // if (dmax < 2) { _beta.clear(); _gamma.clear(); _omega.clear(); return; }

    /* ------------------------------------------------------------
     * 1. cluster volumes  Vol_c = Σ_{v∈c} nodeVolume(v)
     *    !!! No weights
     * ---------------------------------------------------------- */
    PartitionID L = 0;
    for (HypernodeID v : nodes()) // get maximum cluster label
      L = std::max<PartitionID>(L, communityID(v));
    ++L; // clusters are 0-based

    std::vector<double> ClusVol(L, 0.0);
    for (HypernodeID v : nodes())
      // ClusVol[communityID(v)] += static_cast<double>(nodeDegree(v));
      // [mariia's suggestion]:
      ClusVol[communityID(v)] += static_cast<double>(nodeWeightedDegree(v));

    // const double vol_H = initialTotalVertexDegree();
    // [mariia's suggestion]:
    const double vol_H = totalVolume();

    /* ------------------------------------------------------------
     * 2. count edges and cut edges per size k
     * ---------------------------------------------------------- */
    // m_k - sum of edge weights per size k
    std::vector<double> m_k(dmax + 1, 0.01); // small bias avoids log(0)
    // cut_k - sum of cut edge weights per size k
    std::vector<double> cut_k(dmax + 1, 0.0);

    for (HyperedgeID e : edges()) {
      const std::size_t d = static_cast<std::size_t>(edgeSize(e));
      if (d < 2) // ignore single pin nets
        continue;
      const double w = static_cast<double>(edgeWeight(e));

      bool cutting = false;
      PartitionID first_c = communityID(*pins(e).begin());
      for (HypernodeID pin : pins(e))
        if (communityID(pin) != first_c) {
          cutting = true;
          break;
        }

      m_k[d] += w;
      if (cutting)
        cut_k[d] += w;
    }

    /* ------------------------------------------------------------
     * 3. turn counts into ω, β, γ
     * ---------------------------------------------------------- */
    _beta.assign(dmax + 1, 0.0);
    _gamma.assign(dmax + 1, 0.0);
    _omega.assign(dmax + 1, {0.0, 0.0});

    for (std::size_t d = 2; d <= dmax; ++d) {
      // sum of d-th powers of community volumes (no weights!!!)
      double vol_in = 0.0; 
      for (double vc : ClusVol)
        vol_in += std::pow(vc, static_cast<int>(d));
      
      const double vol_out = std::pow(vol_H, static_cast<int>(d)) - vol_in;

      // double omega_in  = (m_k[d] - cut_k[d]) / std::max(vol_in , eps);
      // double omega_out =  cut_k[d]            / std::max(vol_out, eps);

      // omega_in  = std::max(omega_in , eps);
      // omega_out = std::max(omega_out, eps);

      double omega_in = (m_k[d] - cut_k[d]) / vol_in;
      double omega_out = cut_k[d] / vol_out;

      _omega[d] = {omega_in, omega_out};
      // _beta[d] = std::log(omega_in / omega_out);
      _beta[d] = std::log(omega_in) -  std::log(omega_out); // [mariia: this is -\beta_k from (15)]
      if (!std::isfinite(_beta[d])) {
        _beta[d] = (_beta[d] > 0 ? 1e3 : -1e3);
      }
      _gamma[d] = omega_in - omega_out; // [mariia: this is -\beta_k * \gamma_k from (15)]
      LOG << "For edge size d = " << d << ": beta_d = " << _beta[d]
          << ", gamma_d = " << _gamma[d];
    }
  }
  // ═══ AON MOD END ═══════════════════════════════════════════════════════

  // ####################### Fixed Vertex Support #######################

  void addFixedVertexSupport(FixedVertexSupport<StaticHypergraph>&& fixed_vertices) {
    _fixed_vertices = std::move(fixed_vertices);
    _fixed_vertices.setHypergraph(this);
  }

  bool hasFixedVertices() const {
    return _fixed_vertices.hasFixedVertices();
  }

  HypernodeWeight totalFixedVertexWeight() const {
    return _fixed_vertices.totalFixedVertexWeight();
  }

  HypernodeWeight fixedVertexBlockWeight(const PartitionID block) const {
    return _fixed_vertices.fixedVertexBlockWeight(block);
  }

  bool isFixed(const HypernodeID hn) const {
    return _fixed_vertices.isFixed(hn);
  }

  PartitionID fixedVertexBlock(const HypernodeID hn) const {
    return _fixed_vertices.fixedVertexBlock(hn);
  }

  void setMaxFixedVertexBlockWeight(const std::vector<HypernodeWeight> max_block_weights) {
    _fixed_vertices.setMaxBlockWeight(max_block_weights);
  }

  const FixedVertexSupport<StaticHypergraph>& fixedVertexSupport() const {
    return _fixed_vertices;
  }

  FixedVertexSupport<StaticHypergraph> copyOfFixedVertexSupport() const {
    return _fixed_vertices.copy();
  }

  // #################### Flags for PHG (here for now) #####################

  // ! To enable collective sync updates in phg
  void enableCollectiveSyncUpdates() {
    /// [debug] std::cerr << "StaticHypergraph::enableCollectiveSyncUpdates()" << std::endl;
    _enable_collective_sync_update = true;
  }
  bool areCollectiveSyncUpdatesEnabled() const {
    /// [debug] std::cerr << "StaticHypergraph::areCollectiveSyncUpdatesEnabled()" << std::endl;
    return _enable_collective_sync_update;
  }

  // ####################### Single-Pin Nets Removal #######################

  // ! Disable single-pin nets removal before first contraction
  // ! Needed for both Objective::conductance_local ans Objective::conductance_global
  void disableSinglePinNetsRemoval() {
    /// [debug] std::cerr << "StaticHypergraph::disableSinglePinNetsRemoval()" << std::endl;
    _disable_single_pin_nets_removal = true;
  }

  // ! Check if single-pin nets removal is disabled
  bool isSinglePinNetsRemovalDisabled() const {
    return _disable_single_pin_nets_removal;
  }

  // ####################### Contract / Uncontract #######################

  /*!
   * Contracts a given community structure. All vertices with the same label
   * are collapsed into the same vertex. The resulting single-pin (if single-pin 
   * nets removal isn't disabled) and parallel hyperedges are removed from the 
   * contracted graph. The function returns the contracted hypergraph and a 
   * mapping which specifies a mapping from community label (given in 'communities') 
   * to a vertex in the coarse hypergraph.
   *
   * \param communities Community structure that should be contracted
   */
  StaticHypergraph contract(parallel::scalable_vector<HypernodeID>& communities, bool deterministic = false);

  bool registerContraction(const HypernodeID, const HypernodeID) {
    throw UnsupportedOperationException(
      "registerContraction(u, v) is not supported in static hypergraph");
    return false;
  }

  size_t contract(const HypernodeID,
                  const HypernodeWeight max_node_weight = std::numeric_limits<HypernodeWeight>::max()) {
    unused(max_node_weight);
    throw UnsupportedOperationException(
      "contract(v, max_node_weight) is not supported in static hypergraph");
    return 0;
  }

  void uncontract(const Batch&,
                  const UncontractionFunction& case_one_func = NOOP_BATCH_FUNC,
                  const UncontractionFunction& case_two_func = NOOP_BATCH_FUNC) {
    unused(case_one_func);
    unused(case_two_func);
    throw UnsupportedOperationException(
      "uncontract(batch) is not supported in static hypergraph");
  }

  VersionedBatchVector createBatchUncontractionHierarchy(const size_t) {
    throw UnsupportedOperationException(
      "createBatchUncontractionHierarchy(batch_size) is not supported in static hypergraph");
    return { };
  }

  // ####################### Remove / Restore Hyperedges #######################

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    _total_volume -= edgeWeight(he) * edgeSize(he);
    for ( const HypernodeID& pin : pins(he) ) {
      _weighted_degrees[pin] -= edgeWeight(he);
      removeIncidentEdgeFromHypernode(he, pin);
    }
    ++_num_removed_hyperedges;
    disableHyperedge(he);
  }

  /*!
  * Removes a hyperedge from the hypergraph. This includes the removal of he from all
  * of its pins and to disable the hyperedge. Noze, in contrast to removeEdge, this function
  * removes hyperedge from all its pins in parallel.
  *
  * NOTE, this function is not thread-safe and should only be called in a single-threaded
  * setting.
  */
  void removeLargeEdge(const HyperedgeID he) {
    ASSERT(edgeIsEnabled(he), "Hyperedge" << he << "is disabled");
    const size_t incidence_array_start = hyperedge(he).firstEntry();
    const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
    _total_volume -= edgeWeight(he) * edgeSize(he);
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _incidence_array[pos];
      _weighted_degrees[pin] -= edgeWeight(he);
      removeIncidentEdgeFromHypernode(he, pin);
    });
    disableHyperedge(he);
  }

  /*!
   * Restores a large hyperedge previously removed from the hypergraph.
   */
  void restoreLargeEdge(const HyperedgeID& he) {
    ASSERT(!edgeIsEnabled(he), "Hyperedge" << he << "is enabled");
    enableHyperedge(he);
    const size_t incidence_array_start = hyperedge(he).firstEntry();
    const size_t incidence_array_end = hyperedge(he).firstInvalidEntry();
    _total_volume += edgeWeight(he) * edgeSize(he);
    tbb::parallel_for(incidence_array_start, incidence_array_end, [&](const size_t pos) {
      const HypernodeID pin = _incidence_array[pos];
      _weighted_degrees[pin] += edgeWeight(he);
      insertIncidentEdgeToHypernode(he, pin);
    });
  }

  parallel::scalable_vector<ParallelHyperedge> removeSinglePinAndParallelHyperedges() {
    throw UnsupportedOperationException(
      "removeSinglePinAndParallelHyperedges() is not supported in static hypergraph");
    return { };
  }

  void restoreSinglePinAndParallelNets(const parallel::scalable_vector<ParallelHyperedge>&) {
    throw UnsupportedOperationException(
      "restoreSinglePinAndParallelNets(hes_to_restore) is not supported in static hypergraph");
  }

  // ####################### Initialization / Reset Functions #######################

  // ! Reset internal community information
  void copyCommunityIDs(const parallel::scalable_vector<PartitionID>& community_ids) {
    ASSERT(community_ids.size() == UI64(_num_hypernodes));
    doParallelForAllNodes([&](const HypernodeID& hn) {
      _community_ids[hn] = community_ids[hn];
    });
  }

  void setCommunityIDs(ds::Clustering&& communities) {
    ASSERT(communities.size() == initialNumNodes());
    _community_ids = std::move(communities);
  }

  // ! Copy static hypergraph in parallel
  StaticHypergraph copy(parallel_tag_t) const;

  // ! Copy static hypergraph sequential
  StaticHypergraph copy() const;

  // ! Reset internal data structure
  void reset() { }

  // ! Free internal data in parallel
  void freeInternalData() {
    if ( _num_hypernodes > 0 || _num_hyperedges > 0 ) {
      freeTmpContractionBuffer();
    }
    _num_hypernodes = 0;
    _num_hyperedges = 0;
  }

  void freeTmpContractionBuffer() {
    if ( _tmp_contraction_buffer ) {
      delete(_tmp_contraction_buffer);
      _tmp_contraction_buffer = nullptr;
    }
  }

  void memoryConsumption(utils::MemoryTreeNode* parent) const;

    // ! Only for testing
  bool verifyIncidenceArrayAndIncidentNets() {
    throw UnsupportedOperationException(
      "verifyIncidenceArrayAndIncidentNets() not supported in static hypergraph");
    return false;
  }

 private:
  friend class StaticHypergraphFactory;
  template<typename Hypergraph>
  friend class CommunitySupport;
  template <typename Hypergraph,
            typename ConnectivityInformation>
  friend class PartitionedHypergraph;

  // ####################### Hypernode Information #######################

  // ! Accessor for hypernode-related information
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hypernode& hypernode(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return _hypernodes[u];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hypernode& hypernode(const HypernodeID u) {
    return const_cast<Hypernode&>(static_cast<const StaticHypergraph&>(*this).hypernode(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE IteratorRange<IncidentNetsIterator> incident_nets_of(const HypernodeID u,
                                                                                          const size_t pos = 0) const {
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");
    const Hypernode& hn = hypernode(u);
    return IteratorRange<IncidentNetsIterator>(
      _incident_nets.cbegin() + hn.firstEntry() + pos,
      _incident_nets.cbegin() + hn.firstInvalidEntry());
  }

  // ####################### Hyperedge Information #######################

  // ! Accessor for hyperedge-related information
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Hyperedge& hyperedge(const HyperedgeID e) const {
    ASSERT(e <= _num_hyperedges, "Hyperedge" << e << "does not exist");
    return _hyperedges[e];
  }

  // ! To avoid code duplication we implement non-const version in terms of const version
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Hyperedge& hyperedge(const HyperedgeID e) {
    return const_cast<Hyperedge&>(static_cast<const StaticHypergraph&>(*this).hyperedge(e));
  }

  // ####################### Remove / Restore Hyperedges #######################

  // ! Removes hyperedge e from the incident nets of vertex hn
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void removeIncidentEdgeFromHypernode(const HyperedgeID e,
                                                                          const HypernodeID u) {
    using std::swap;
    ASSERT(!hypernode(u).isDisabled(), "Hypernode" << u << "is disabled");

    Hypernode& hn = hypernode(u);
    size_t incident_nets_pos = hn.firstEntry();
    for ( ; incident_nets_pos < hn.firstInvalidEntry(); ++incident_nets_pos ) {
      if ( _incident_nets[incident_nets_pos] == e ) {
        break;
      }
    }
    ASSERT(incident_nets_pos < hn.firstInvalidEntry());
    swap(_incident_nets[incident_nets_pos], _incident_nets[hn.firstInvalidEntry() - 1]);
    hn.setSize(hn.size() - 1);
  }

  // ! Inserts hyperedge he to incident nets array of vertex hn
  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void insertIncidentEdgeToHypernode(const HyperedgeID e,
                                                                        const HypernodeID u) {
    using std::swap;
    Hypernode& hn = hypernode(u);
    ASSERT(!hn.isDisabled(), "Hypernode" << u << "is disabled");
    HEAVY_REFINEMENT_ASSERT(std::count(_incident_nets.cbegin() + hn.firstEntry(),
                                       _incident_nets.cbegin() + hn.firstInvalidEntry(), e) == 0,
                        "HN" << u << "is already connected to HE" << e);
    const size_t incident_nets_start = hn.firstInvalidEntry();
    const size_t incident_nets_end = hypernode(u + 1).firstEntry();
    size_t incident_nets_pos = incident_nets_start;
    for ( ; incident_nets_pos < incident_nets_end; ++incident_nets_pos ) {
      if ( _incident_nets[incident_nets_pos] == e ) {
        break;
      }
    }
    ASSERT(incident_nets_pos < incident_nets_end);
    swap(_incident_nets[incident_nets_start], _incident_nets[incident_nets_pos]);
    hn.setSize(hn.size() + 1);
  }

  // ! Allocate the temporary contraction buffer
  void allocateTmpContractionBuffer() {
    if ( !_tmp_contraction_buffer ) {
      _tmp_contraction_buffer = new TmpContractionBuffer(
        _num_hypernodes, _num_hyperedges, _num_pins);
    }
  }

  // ! Number of hypernodes
  HypernodeID _num_hypernodes;
  // ! Number of removed hypernodes
  HypernodeID _num_removed_hypernodes;
  // ! Number of removed degree zero hypernodes
  HypernodeWeight _removed_degree_zero_hn_weight;
  // ! Number of hyperedges
  HyperedgeID _num_hyperedges;
  // ! Number of removed hyperedges
  HyperedgeID _num_removed_hyperedges;
  // ! Maximum size of a hyperedge
  HypernodeID _max_edge_size;
  // ! Maximum size of a hyperedge at the moment of the last snapshot
  HypernodeID _original_max_edge_size;
  // ! Number of pins
  HypernodeID _num_pins;
  // ! Total degree of all vertices
  HypernodeID _total_degree;
  // ! Total weight of hypergraph
  HypernodeWeight _total_weight;
  // ! Total volume of hypergraph
  HypergraphVolume _total_volume;
  // ! Original total volume of hypergraph
  HypergraphVolume _original_total_volume;

  // ! Hypernodes
  Array<Hypernode> _hypernodes;
  // ! Pins of hyperedges
  IncidentNets _incident_nets;
  // ! Hyperedges
  Array<Hyperedge> _hyperedges;
  // ! Incident nets of hypernodes
  IncidenceArray _incidence_array;
  // ! Weighted degrees of hypernodes
  Array<HypergraphVolume> _weighted_degrees;
  // ! Original Weighted degrees of hypernodes
  // ! (are lost during contractions, removal of single-pin nets)
  Array<HypergraphVolume> _original_weighted_degrees;

  // ! Communities
  ds::Clustering _community_ids;

  // ! Fixed Vertex Support
  FixedVertexSupport<StaticHypergraph> _fixed_vertices;

  // ! Data that is reused throughout the multilevel hierarchy
  // ! to contract the hypergraph and to prevent expensive allocations
  TmpContractionBuffer* _tmp_contraction_buffer;

  // ! Option for enabling sync_updates in phg
  bool _enable_collective_sync_update = false;

  // ! Option for disabling the removal of single-pin nets
  bool _disable_single_pin_nets_removal = false;

  
  // AON HyperModularity Clustering Coefficients
  vec<double> _beta;                 ///< -β_k
  vec<double> _gamma;                ///< -β_k * γ_k
  vec<std::array<double, 2>> _omega; ///< {ω_k0, ω_k1} (ω_in, ω_out)
};

} // namespace ds
} // namespace mt_kahypar
