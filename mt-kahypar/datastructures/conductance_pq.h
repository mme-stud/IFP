
#pragma once

#include <atomic>
#include <mutex>

#include <tbb/parallel_invoke.h>

#include "kahypar-resources/meta/mandatory.h"

#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/priority_queue.h"

namespace mt_kahypar {
namespace ds {

template <typename Numerator, typename Denominator = Numerator>
class NonnegativeFraction {
private:
  Numerator numerator;
  Denominator denominator;
public:
  NonnegativeFraction() :
    numerator(0),
    denominator(1) { }

  NonnegativeFraction(const Numerator& n, const Denominator& d) :
    numerator(n),
    denominator(d) { 
      ASSERT(d > 0);
      ASSERT(n > 0);
    }

  bool operator< (const NonnegativeFraction& other) const {
    size_t lhs = static_cast<size_t>(numerator) * static_cast<size_t>(other.denominator);
    size_t rhs = static_cast<size_t>(other.numerator) * static_cast<size_t>(denominator);
    return lhs < rhs;
  }

  bool operator== (const NonnegativeFraction& other) const {
    size_t lhs = static_cast<size_t>(numerator) * static_cast<size_t>(other.denominator);
    size_t rhs = static_cast<size_t>(other.numerator) * static_cast<size_t>(denominator);
    return lhs == rhs;
  }

  bool operator> (const NonnegativeFraction& other) const {
    size_t lhs = static_cast<size_t>(numerator) * static_cast<size_t>(other.denominator);
    size_t rhs = static_cast<size_t>(other.numerator) * static_cast<size_t>(denominator);
    return lhs > rhs;
  }

  // operator<< 
  friend std::ostream& operator<< (std::ostream& s, const NonnegativeFraction& f) {
    return s << f.numerator << " / " << f.denominator;
  }

  // Numerator must be non-negative 
  void setNumerator(const Numerator& n) {
    ASSERT(n >= 0);
    numerator = n;
  }
  // ! Denominator must be greater than 0
  void setDenominator(const Denominator& d) {
    ASSERT(d > 0);
    denominator = d;
  }
  Numerator getNumerator() const {
    return numerator;
  }
  Denominator getDenominator() const {
    return denominator;
  }

  double_t value() const {
    return static_cast<double_t>(numerator) / static_cast<double_t>(denominator);
  }
};

using ConductanceFraction = NonnegativeFraction<HyperedgeWeight>;


/**
 * @brief Priority Queue for partitions based on their conductance.
 * 
 * All write operations are synchronized.
 * All read operations aren't synchronized per default (use the synchronized flag to synchronize).
 */
template <typename PartitionedHypergraph = Mandatory>
class ConductancePriorityQueue : 
      protected ExclusiveHandleHeap<MaxHeap<ConductanceFraction, PartitionID>> {
private:
  using SuperPQ = ExclusiveHandleHeap<MaxHeap<ConductanceFraction, PartitionID>>;
public:
  ConductancePriorityQueue() :
    SuperPQ(0),
    _total_volume(-1),
    _size(0),
    _complement_val_bits(),
    _initialized(false)
    { }
  
  // ! Initializes the priority queue with the partitions of the hypergraph
  void initialize(const PartitionedHypergraph& hg, bool synchronized = false) {
    lock(synchronized);
    ASSERT(!_initialized);
    _initialized = true;
    _total_volume = getHGTotalVolume(hg);
    _size = hg.k();
    _complement_val_bits.resize(_size);
    SuperPQ::clear();
    SuperPQ::resize(_size);
    SuperPQ::heap.resize(_size);
    tbb::parallel_for(PartitionID(0), _size, [&](const PartitionID& p) {
      HyperedgeWeight cut_weight = getHGPartCutWeight(hg, p);
      HyperedgeWeight volume = getHGPartVolume(hg, p);
      _complement_val_bits[p] = (volume > _total_volume - volume);
      ConductanceFraction f(cut_weight, std::min(volume, _total_volume - volume));
      SuperPQ::heap[p].id = p; 
      SuperPQ::heap[p].key = f;
      SuperPQ::positions[p] = p;
    });
    buildHeap();
    unlock(synchronized);
  }

  // ! Returns true if the priority queue is initialized
  bool initialized() const {
    return _initialized;
  }

  // ! Reset the priority queue to the uninitialized state
  void reset(bool synchronized = false) {
    lock(synchronized);
    SuperPQ::clear();
    _total_volume = -1;
    _size = 0;
    _complement_val_bits.clear();
    _initialized = false;
    unlock(synchronized);
  }

  // ! Returns an approximate memory consumption of the conductance priority queue in bytes
  size_t memoryConsumption() const {
    return SuperPQ::memoryConsumption() + _complement_val_bits.size() * sizeof(bool);
  }

  // ! Updates the priority queue after global changes in partition
  void globalUpdate(const PartitionedHypergraph& hg, bool synchronized = false) {
    lock(synchronized);
    ASSERT(_initialized && _size == hg.k());
    if (_uses_original_stats) {
      ASSERT(_total_volume == hg.originalTotalVolume());
    }
    _total_volume = getHGTotalVolume(hg);
    tbb::parallel_for(PartitionID(0), _size, [&](const PartitionID& p) {
      HyperedgeWeight cut_weight = getHGPartCutWeight(hg, p);
      HyperedgeWeight volume = getHGPartVolume(hg, p);
      _complement_val_bits[p] = (volume > _total_volume - volume);
      ConductanceFraction f(cut_weight, std::min(volume, _total_volume - volume));
      SuperPQ::heap[SuperPQ::positions[p]].key = f;
    });
    buildHeap();
    unlock(synchronized);
  }

  // ! Checks if the priority queue is correct with respect to the hypergraph: const version
  bool check(const PartitionedHypergraph& hg) const {
    ASSERT(_initialized && _size == hg.k());
    bool correct = true;
    if (_total_volume != getHGTotalVolume(hg)) {
      correct = false;
      LOG << "Total volume in ConductancePriorityQueue is" << _total_volume << ", but should be" << getHGTotalVolume(hg);
    }
    for (PartitionID p = 0; p < _size; ++p) {
      ConductanceFraction f = SuperPQ::getKey(p);
      HyperedgeWeight cut_weight = f.getNumerator();
      HyperedgeWeight volume = f.getDenominator();
      if (_complement_val_bits[p]) {
        volume = _total_volume - volume;
      }
      if (volume != getHGPartVolume(hg, p)) {
        correct = false;
        LOG << "Volume of partition in ConductancePriorityQueue" << p << "is" << volume << ", but should be" << getHGPartVolume(hg, p);
      }
      if (cut_weight != getHGPartCutWeight(hg, p)) {
        correct = false;
        LOG << "Cut weight of partition in ConductancePriorityQueue" << p << "is" << cut_weight << ", but should be" << getHGPartCutWeight(hg, p);
      }
    }
    correct = correct && SuperPQ::isHeap() && SuperPQ::positionsMatch();
    return correct;
  }

  // ! Checks if the priority queue is correct with respect to the hypergraph: synchronizable version
  bool checkSync(const PartitionedHypergraph& hg, bool synchronized = false) {
    lock(synchronized);
    bool correct = check(hg);
    unlock(synchronized);
    return correct;
  }

  // ################# Priority Queue Operations #################

  // ! Adjusts the cut weight and the volume of a partition
  // ! changes pq => uses a lock  
  void adjustKey(const PartitionID& p, const HyperedgeWeight& cut_weight, const HyperedgeWeight& volume, bool synchronized = true) {
    ASSERT(_total_volume >= volume && volume >= 0);
    _complement_val_bits[p] = (volume > _total_volume - volume);
    ConductanceFraction f(cut_weight, std::min(volume, _total_volume - volume));
    lock(synchronized);
    SuperPQ::adjustKey(p, f);
    unlock(synchronized);
  }
  
  // ! Updates PQ after total volume of the hypergraph has changed
  // ! changes pq => uses a lock
  void updateTotalVolume(const HyperedgeWeight& new_total_volume, bool synchronized = true) {
    ASSERT(!_use_original_stats);
    lock(synchronized);
    for (PartitionID p = 0; p < _size; ++p) {
      ConductanceFraction f = SuperPQ::getKey(p);
      HyperedgeWeight volume = f.getDenominator();
      if (_complement_val_bits[p]) {
        volume = _total_volume - volume;
      }
      f.setDenominator(std::min(volume, new_total_volume - volume));
    }
    buildHeap();
    _total_volume = new_total_volume;
    unlock(synchronized);
  }

  // ! Get the partition with the highest conductance: const version
  PartitionID top() const {
    PartitionID p = SuperPQ::top();
    return p;
  }

  // ! Get the partition with the highest conductance: synchronizable version
  PartitionID topSync(bool synchronized = false) {
    lock(synchronized);
    PartitionID p = top();
    unlock(synchronized);
    return p;
  }

  // ! Get the partition with the second highest conductance: const version
  PartitionID secondTop() const {
    ASSERT(SuperPQ::size() > 1);
    PartitionID f = SuperPQ::heap[1].id;
    // ConductancePriorityQueue is a MaxHeap => binary tree
    if (size() > 2 && SuperPQ::heap[1].key < SuperPQ::heap[2].key) {
      f = SuperPQ::heap[2].id;
    }
    return f;
  }
  
  // ! Get the partition with the second highest conductance: synchronizable version
  PartitionID secondTopSync(bool synchronized = false) {
    lock(synchronized);
    PartitionID f = secondTop();
    unlock(synchronized);
    return f;
  }

  // ! Get the top three partitions (unsorted): constant version
  // ! (Works only for a binary heap)
  vec<PartitionID> topThree() const {
    vec<PartitionID> top_three(3, kInvalidPartition);
    for (size_t i = 0; i < 3 && i < _size; ++i) {
      top_three[i] = SuperPQ::heap[i].id;
    }
    return top_three;
  }

  // ! Get the top three partitions (unsorted): synchronizable version
  // ! (Works only for a binary heap)
  vec<PartitionID> topThreeSync(bool synchronized = false) {
    lock(synchronized);
    vec<PartitionID> top_three = topThree();
    unlock(synchronized);
    return top_three;
  }
  
  bool empty() const {
    bool e = SuperPQ::empty();
    return e;
  }

  size_t size() const {
    ASSERT(static_cast<PosT>(_size) == SuperPQ::size());
    return _size;
  }

  // ################## USAGE OF ORIGINAL PHG STATS #################

  // ! Uses the original stats of the hypergraph
  bool usesOriginalStats() const {
    return _uses_original_stats;
  }
  
  // ! Makes ConductancePriorityQueue use current _total_volume and _part_volumes
  // ! instead of the original (inherited from old versions of hg) ones
  // ! To be used before initialization
  void disableUsageOfOriginalHGStats() {
    ASSERT(!_initialized, "ConductancePriorityQueue is already initialized");
    _uses_original_stats = false;
  }

  // ! Makes ConductancePriorityQueue use original _total_volume and _part_volumes
  // ! instead of the current ones
  // ! To be used before initialization
  void enableUsageOfOriginalHGStats() {
    ASSERT(!_initialized, "ConductancePriorityQueue is already initialized");
    _uses_original_stats = true;
  }

private:
  // ! Builds the heap in O(_size) time
  // ! no built in lock
  void buildHeap() {
    ASSERT(static_cast<PosT>(_size) == SuperPQ::size());
    if (SuperPQ::isHeap()) return;
    for (PartitionID p = _size - 1; p >= 0; --p) {
      SuperPQ::siftDown(p);
    }
    ASSERT(SuperPQ::isHeap() && SuperPQ::positionsMatch());
  }

  // ################### COMMUNICATION WITH THE HG ######################
  // ! Get needed kind of total volume
  // ! (original or current)
  HyperedgeWeight getHGTotalVolume(const PartitionedHypergraph& hg) {
    if (_uses_original_stats) {
      return hg.originalTotalVolume();
    } else {
      return hg.totalVolume();
    }
  }
  
  // ! Get needed kind of part volume
  // ! (original or current)
  HyperedgeWeight getHGPartVolume(const PartitionedHypergraph& hg, const PartitionID p) {
    if (_uses_original_stats) {
      return hg.partOriginalVolume(p);
    } else {
      return hg.partVolume(p);
    }
  }

  // ! Get needed kind of cut weight
  // ! (only one kind of cut weight for now)
  HyperedgeWeight getHGPartCutWeight(const PartitionedHypergraph& hg, const PartitionID p) {
    return hg.partCutWeight(p);
  }

  // #################### POTENTIALLY USELESS PART ######################
  
  // ! insert a partition with its cut weight and volume
  // ! changes pq => uses a lock
  // ! not sure if this is useful
  void insert(const PartitionID& p, const HyperedgeWeight& cut_weight, const HyperedgeWeight& volume, bool synchronized = true) {
    ConductanceFraction f(cut_weight, std::min(volume, _total_volume - volume));
    lock(synchronized);
    SuperPQ::insert(p, f);
    _size++;
    unlock(synchronized);
  }

  // ! changes pq => uses a lock
  // ! not sure if this is useful
  void remove(const PartitionID& p, bool synchronized = true) {
    lock(synchronized);
    SuperPQ::remove(p);
    _size--;
    unlock(synchronized);
  }


  // ! not sure if this is useful
  void deleteTop(bool synchronized = true) {
    lock(synchronized);
    SuperPQ::deleteTop();
    _size--;
    unlock(synchronized);
  }

  // ! Get the conductance fraction of the partition with the highest conductance
  ConductanceFraction topFractionSync(bool synchronized = false) {
    lock(synchronized);
    ConductanceFraction f = SuperPQ::topKey();
    unlock(synchronized);
    return f;
  }
  // ! Get the conductance fraction of the partition with the second highest conductance
  ConductanceFraction secondTopFractionSync(bool synchronized = false) {
    ASSERT(SuperPQ::size() > 1);
    lock(synchronized);
    ConductanceFraction f =  SuperPQ::heap[1].key;
    // ConductancePriorityQueue is a MaxHeap => binary tree
    if (SuperPQ::size() > 2) {
      f = std::max(f, SuperPQ::heap[2].key);
    }
    unlock(synchronized);
    return f;
  }

  double_t topConductanceSync(bool synchronized = false) {
    ConductanceFraction f = topFractionSync(synchronized);
    return f.value();
  }
  double_t secondTopCondunctanceSync(bool synchronized = false) {
    ConductanceFraction f = secondTopFractionSync(synchronized);
    return f.value();
  }

  // ! Get the conductance fraction of a partition
  // (better use PartitionedHypergraph getPartCutWeight and getPartVolume)
  ConductanceFraction getFractionSync(const PartitionID p, bool synchronized = false) {
    lock(synchronized);
    ConductanceFraction f = SuperPQ::getKey(p);
    unlock(synchronized);
    return f;
  }
  // ! Get the conductance of a partition
  // (better use PartitionedHypergraph getPartCutWeight and getPartVolume)
  double_t getConductanceSync(const PartitionID p, bool synchronized = false) {
    lock(synchronized);
    ConductanceFraction f = SuperPQ::getKey(p);
    unlock(synchronized);
    return f.value();
  }

  // ################# MUTEX OPERATIONS FOR CHANGING PQ #################

  void lock(bool synchronized) {
    if (synchronized) _pq_lock.lock();
  }

  void unlock(bool synchronized) {
    if (synchronized) _pq_lock.unlock();
  }

  // ################# MEMBER VARIABLES #################
  SpinLock _pq_lock;
  HyperedgeWeight _total_volume;
  PartitionID _size;
  vec<bool> _complement_val_bits;
  bool _uses_original_stats = true;
  bool _initialized;
};


}  // namespace ds
}  // namespace mt_kahypar