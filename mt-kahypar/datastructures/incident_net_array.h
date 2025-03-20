/*******************************************************************************
 * MIT License
 *
 * This file is part of Mt-KaHyPar.
 *
 * Copyright (C) 2020 Lars Gottesbüren <lars.gottesbueren@kit.edu>
 * Copyright (C) 2020 Tobias Heuer <tobias.heuer@kit.edu>
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

#include <cstddef>

#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_invoke.h>
#include <tbb/parallel_scan.h>

#include "kahypar-resources/datastructure/fast_reset_flag_array.h"

#include "mt-kahypar/macros.h"
#include "mt-kahypar/datastructures/hypergraph_common.h"
#include "mt-kahypar/datastructures/array.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/parallel/stl/scalable_unique_ptr.h"
#include "mt-kahypar/parallel/atomic_wrapper.h"
#include "mt-kahypar/parallel/parallel_prefix_sum.h"
#include "mt-kahypar/utils/range.h"

// #include "mt-kahypar/datastructures/dynamic_hypergraph.h" // to calculate weighted_degrees
class DynamicHypergraph;

namespace mt_kahypar {
namespace ds {

// forward declaration
class IncidentNetArray;

// Iterator over the incident nets of a vertex u
class IncidentNetIterator {
  public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = HyperedgeID;
  using reference = HyperedgeID&;
  using pointer = const HyperedgeID*;
  using difference_type = std::ptrdiff_t;

  IncidentNetIterator(const HypernodeID u,
                      const IncidentNetArray* incident_net_array,
                      const size_t pos,
                      const bool end);

  HyperedgeID operator* () const;

  IncidentNetIterator & operator++ ();

  IncidentNetIterator operator++ (int) {
    IncidentNetIterator copy = *this;
    operator++ ();
    return copy;
  }

  bool operator!= (const IncidentNetIterator& rhs);

  bool operator== (const IncidentNetIterator& rhs);

  private:
  void next_iterator();

  HypernodeID _u;
  HypernodeID _current_u;
  HypernodeID _current_size;
  size_t _current_pos;
  const IncidentNetArray* _incident_net_array;
  bool _end;
};

// ! Class allows in-place contraction and uncontraction of the incident net array
class IncidentNetArray {

  using HyperedgeVector = parallel::scalable_vector<parallel::scalable_vector<HypernodeID>>;
  using ThreadLocalCounter = tbb::enumerable_thread_specific<parallel::scalable_vector<size_t>>;
  using AtomicCounter = parallel::scalable_vector<parallel::IntegralAtomicWrapper<size_t>>;

  using AcquireLockFunc = std::function<void (const HypernodeID)>;
  using ReleaseLockFunc = std::function<void (const HypernodeID)>;
  using CaseOneFunc = std::function<void (const HyperedgeID)>;
  using CaseTwoFunc = std::function<void (const HyperedgeID)>;
  #define NOOP_LOCK_FUNC [] (const HypernodeID) { }

  static_assert(sizeof(char) == 1);

  // Represents one incident net of a vertex.
  // A incident net is associated with a version number. Incident nets
  // with a version number greater or equal than the version number in
  // header (see Header -> current_version) are active.
  struct Entry {
    HyperedgeID e;
    HypernodeID version;
  };

  // Header of the incident net list of a vertex. The incident net lists
  // contracted into one vertex are concatenated in a double linked list.
  struct Header {
    explicit Header(const HypernodeID u) :
      prev(u),
      next(u),
      it_prev(u),
      it_next(u),
      tail(u),
      size(0),
      degree(0),
      weighted_degree(0),
      current_version(0),
      is_head(true) { }

    // ! Previous incident net list
    HypernodeID prev;
    // ! Next incident net list
    HypernodeID next;
    // ! Previous non-empty incident net list
    HypernodeID it_prev;
    // ! Next non-empty incident net list
    HypernodeID it_next;
    // ! If we append a vertex v to the incident net list of a vertex u, we store
    // ! the previous tail of vertex v, such that we can restore the list of v
    // ! during uncontraction
    HypernodeID tail;
    // ! All incident nets between [0,size) are active
    HypernodeID size;
    // ! Degree of the vertex
    HypernodeID degree;
    // ! Weighted degree of the vertex (calculated only if _hypergraph_ref for the
    // ! incident net array is not nullptr)
    HyperedgeWeight weighted_degree;
    // ! Current version of the incident net list
    HypernodeID current_version;
    // ! True, if the vertex is the head of a incident net list
    bool is_head;
  };

 public:
  using const_iterator = IncidentNetIterator;

  IncidentNetArray() :
    _num_hypernodes(0),
    _size_in_bytes(0),
    _index_array(),
    _incident_net_array(nullptr),
    _hypergraph_ptr(nullptr) { 
      /// [debug] std::cerr << "IncidentNetArray() !!!" << std::endl;
    }

  IncidentNetArray(const HypernodeID num_hypernodes,
                   const HyperedgeVector& edge_vector,
                   DynamicHypergraph* hypergraph_ptr = nullptr,
                   const HyperedgeWeight* hyperedge_weight_ptr = nullptr) :
    _num_hypernodes(num_hypernodes),
    _size_in_bytes(0),
    _index_array(),
    _incident_net_array(nullptr),
    _hypergraph_ptr(hypergraph_ptr) {
    construct(edge_vector, hyperedge_weight_ptr);
  }

  // ! Degree of the vertex
  HypernodeID nodeDegree(const HypernodeID u) const {
    /// [debug] std::cerr << "nodeDegree(u)" << std::endl;
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return header(u)->degree;
  }

  // ! Weighted degree of the vertex
  HyperedgeWeight nodeWeightedDegree(const HypernodeID u) const {
    /// [debug] std::cerr << "nodeWeightedDegree(u)" << std::endl;
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return header(u)->weighted_degree;
  }

  // ! Decrease weighted degree of the vertex (for updating weighted degrees)
  void decreaseNodeWeightedDegree(const HypernodeID u, HyperedgeWeight w) const {
    /// [debug] std::cerr << "decreaseNodeWeightedDegree(u, w)" << std::endl;
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    // header(u)->weighted_degree -= w;
    Header* headerU = reinterpret_cast<Header*>(_incident_net_array.get() + _index_array[u]);
    headerU->weighted_degree -= w;
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetIterator> incidentEdges(const HypernodeID u) const {
    /// [debug] std::cerr << "incidentEdges(u)" << std::endl;
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentNetIterator>(
      IncidentNetIterator(u, this, UL(0), false),
      IncidentNetIterator(u, this, UL(0), true));
  }

  // ! Returns a range to loop over the incident nets of hypernode u.
  IteratorRange<IncidentNetIterator> incidentEdges(const HypernodeID u,
                                                   const size_t pos) const {
    /// [debug] std::cerr << "incidentEdges(u, pos)" << std::endl;
    ASSERT(u < _num_hypernodes, "Hypernode" << u << "does not exist");
    return IteratorRange<IncidentNetIterator>(
      IncidentNetIterator(u, this, pos, false),
      IncidentNetIterator(u, this, UL(0), true));
  }

  // ! Contracts two incident list of u and v, whereby u is the representative and
  // ! v the contraction partner of the contraction. The contraction involves to remove
  // ! all incident nets shared between u and v from the incident net list of v and append
  // ! the list of v to u.
  void contract(const HypernodeID u,
                const HypernodeID v,
                const kahypar::ds::FastResetFlagArray<>& shared_hes_of_u_and_v,
                const AcquireLockFunc& acquire_lock = NOOP_LOCK_FUNC,
                const ReleaseLockFunc& release_lock = NOOP_LOCK_FUNC);

  // ! Uncontract two previously contracted vertices u and v.
  // ! Uncontraction involves to decrement the version number of all incident lists contained
  // ! in v and restore all incident nets with a version number equal to the new version.
  // ! Note, uncontraction must be done in relative contraction order
  void uncontract(const HypernodeID u,
                  const HypernodeID v,
                  const AcquireLockFunc& acquire_lock = NOOP_LOCK_FUNC,
                  const ReleaseLockFunc& release_lock = NOOP_LOCK_FUNC);

  // ! Uncontract two previously contracted vertices u and v.
  // ! Uncontraction involves to decrement the version number of all incident lists contained
  // ! in v and restore all incident nets with a version number equal to the new version.
  // ! Additionally it calls case_one_func for a hyperedge he, if u and v were previously both
  // ! adjacent to he and case_two_func if only v was previously adjacent to he.
  // ! Note, uncontraction must be done in relative contraction order
  void uncontract(const HypernodeID u,
                  const HypernodeID v,
                  const CaseOneFunc& case_one_func,
                  const CaseTwoFunc& case_two_func,
                  const AcquireLockFunc& acquire_lock,
                  const ReleaseLockFunc& release_lock);

  // ! Removes all incidents nets of u flagged in hes_to_remove.
  void removeIncidentNets(const HypernodeID u,
                          const kahypar::ds::FastResetFlagArray<>& hes_to_remove,
                          bool update_weighted_degrees = true);

  // ! Restores all previously removed incident nets
  // ! Note, function must be called in reverse order of calls to
  // ! removeIncidentNets(...) and all uncontraction that happens
  // ! between two consecutive calls to removeIncidentNets(...) must
  // ! be processed.
  void restoreIncidentNets(const HypernodeID u,
                          bool update_weighted_degrees = true);

  IncidentNetArray copy(parallel_tag_t, DynamicHypergraph* _hypergraph_ptr = nullptr) const;

  IncidentNetArray copy(DynamicHypergraph* _hypergraph_ptr = nullptr) const;

  void reset();

  size_t size_in_bytes() const {
    return _size_in_bytes + sizeof(size_t) * _index_array.size();
  }

 private:
  friend class IncidentNetIterator;

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Header* header(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Header*>(_incident_net_array.get() + _index_array[u]);
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Header* header(const HypernodeID u) {
    return const_cast<Header*>(static_cast<const IncidentNetArray&>(*this).header(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Entry* firstEntry(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Entry*>(_incident_net_array.get() + _index_array[u] + sizeof(Header));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Entry* firstEntry(const HypernodeID u) {
    return const_cast<Entry*>(static_cast<const IncidentNetArray&>(*this).firstEntry(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE const Entry* lastEntry(const HypernodeID u) const {
    ASSERT(u <= _num_hypernodes, "Hypernode" << u << "does not exist");
    return reinterpret_cast<const Entry*>(_incident_net_array.get() +
      _index_array[u] + sizeof(Header) + header(u)->size * sizeof(Entry));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE Entry* lastEntry(const HypernodeID u) {
    return const_cast<Entry*>(static_cast<const IncidentNetArray&>(*this).lastEntry(u));
  }

  MT_KAHYPAR_ATTRIBUTE_ALWAYS_INLINE void swap(Entry* lhs, Entry* rhs) {
    Entry tmp_lhs = *lhs;
    *lhs = *rhs;
    *rhs = tmp_lhs;
  }

  // ! Restores all previously removed incident nets
  // ! Note, function must be called in reverse order of calls to
  // ! removeIncidentNets(...) and all uncontraction that happens
  // ! between two consecutive calls to removeIncidentNets(...) must
  // ! be processed.
  void restoreIncidentNets(const HypernodeID u,
                           const CaseOneFunc& case_one_func,
                           const CaseTwoFunc& case_two_func,
                           bool update_weighted_degrees = true);

  void append(const HypernodeID u, const HypernodeID v);

  void splice(const HypernodeID u, const HypernodeID v);

  void removeEmptyIncidentNetList(const HypernodeID u);

  void construct(const HyperedgeVector& edge_vector, const HyperedgeWeight* hyperedge_weight_ptr = nullptr);

  bool verifyIteratorPointers(const HypernodeID u) const;

  HypernodeID _num_hypernodes;
  size_t _size_in_bytes;
  Array<size_t> _index_array;
  parallel::tbb_unique_ptr<char> _incident_net_array;
  DynamicHypergraph* _hypergraph_ptr; // to calculate weighted_degrees
};

}  // namespace ds
}  // namespace mt_kahypar