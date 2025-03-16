
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
class Fraction {
private:
  Numerator numerator;
  Denominator denominator;
public:
  Fraction() :
    numerator(0),
    denominator(1) { }

  Fraction(const Numerator& n, const Denominator& d) :
    numerator(n),
    denominator(d) { }

  Fraction(const Fraction& other) :
    numerator(other.numerator),
    denominator(other.denominator) { }

  bool operator< (const Fraction& other) const {
    size_t lhs = static_cast<size_t>(numerator) * static_cast<size_t>(other.denominator);
    size_t rhs = static_cast<size_t>(other.numerator) * static_cast<size_t>(denominator);
    return lhs < rhs;
  }

  void setNumerator(const Numerator& n) {
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

using ConductanceFraction = Fraction<HyperedgeWeight>;


/**
 * @brief Priority Queue for partitions based on their conductance
 * All write operations are synchronized
 * All read operations aren't synchronized per default (use the synchronized flag to synchronize)
 */
template <typename PartitionedHypergraph = Mandatory>
class ConductancePriorityQueue : 
      protected ExclusiveHandleHeap<MaxHeap<PartitionID, ConductanceFraction>> {
private:
  using SuperPQ = ExclusiveHandleHeap<MaxHeap<PartitionID, ConductanceFraction>>;
public:
  ConductancePriorityQueue() :
    SuperPQ(0),
    _total_volume(-1),
    _size(0),
    _complement_val_bits(),
    _initialized(false)
    { }
  
  ConductancePriorityQueue(const ConductancePriorityQueue& other) :
    SuperPQ(other),
    _total_volume(other._total_volume),
    _size(other._size),
    _complement_val_bits(other._complement_val_bits),
    _initialized(other._initialized)
    { }
  
  // ! Initializes the priority queue with the partitions of the hypergraph
  void initialize(const PartitionedHypergraph& hg) {
    ASSERT(!_initialized);
    _initialized = true;
    _total_volume = hg.totalVolume();
    _size = hg.k();
    _complement_val_bits.resize(_size);
    SuperPQ::clear();
    SuperPQ::resize(_size);
    SuperPQ::heap.resize(_size);
    tbb::parallel_for(PartitionID(0), _size, [&](const PartitionID& p) {
      HyperedgeWeight cut_weight = hg.partCutWeight(p);
      HyperedgeWeight volume = hg.partVolume(p);
      _complement_val_bits[p] = (volume > _total_volume - volume);
      ConductanceFraction f(cut_weight, std::min(volume, _total_volume - volume));
      SuperPQ::heap[p] = {p, f};
      SuperPQ::positions[p] = p;
    });
    buildHeap();
  }

  // ! Returns true if the priority queue is initialized
  bool initialized() const {
    return _initialized;
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
    lock(synchronized);
    for (PartitionID p = 0; p < _size; ++p) {
      ConductanceFraction f = SuperPQ::getKey(p);
      HyperedgeWeight volume = f.getDenominator();
      if (complement_val_bits[p]) {
        volume = _total_volume - volume;
      }
      f.setDenominator(std::min(volume, new_total_volume - volume));
    }
    buildHeap();
    _total_volume = new_total_volume;
    unlock(synchronized);
  }

  // ! Get the partition with the highest conductance
  PartitionID top(bool syncronized = false) const {
    lock(syncronized);
    PartitionID p = SuperPQ::top();
    unlock(syncronized);
    return p;
  }
  // ! Get the partition with the second highest conductance
  PartitionID secondTop(bool syncronized = false) const {
    ASSERT(SuperPQ::size() > 1);
    lock(syncronized);
    PartitionID f =  SuperPQ::heap[1].id;
    // ConductancePriorityQueue is a MaxHeap => binary tree
    if (size() > 2 && SuperPQ::heap[1].value < SuperPQ::heap[2].value) {
      f = SuperPQ::heap[2].id;
    }
    unlock(syncronized);
    return f;
  }

  // ! Get the top three partitions (unsorted)
  // ! (Works only for a binary heap)
  vec<PartitionID> topThree(bool synchronized = false) const {
    lock(synchronized);
    vec<PartitionID> topThree(3, kInvalidPartition);
    for (size_t i = 0; i < 3 && i < _size; ++i) {
      topThree[i] = SuperPQ::heap[i].id;
    }
    unlock(synchronized);
    return topThree;
  }
  
  bool empty() const {
    bool e = SuperPQ::empty();
    return e;
  }

  size_t size() const {
    ASSERT(_k == SuperPQ::size());
    return _k;
  }

private:
  // ! Builds the heap in O(_size) time
  // ! no built in lock
  void buildHeap() {
    ASSERT(_size == SuperPQ::size());
    if SuperPQ::isHeap() return;
    for (PartitionID p = _size - 1; p >= 0; --p) {
      SuperPQ::siftDown(p);
    }
    ASSERT(SuperPQ::isHeap() && SuperPQ::positionsMatch());
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
  void clear(bool synchronized = true) {
    lock(synchronized);
    SuperPQ::clear();
    _total_volume = -1;
    _size = 0;
    _complement_val_bits.clear();
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
  ConductanceFraction topFraction(bool syncronized = false) const {
    lock(syncronized);
    ConductanceFraction f = SuperPQ::topKey();
    unlock(syncronized);
    return f;
  }
  // ! Get the conductance fraction of the partition with the second highest conductance
  ConductanceFraction secondTopFraction(bool syncronized = false) const {
    ASSERT(SuperPQ::size() > 1);
    lock(syncronized);
    ConductanceFraction f =  SuperPQ::heap[1].value;
    // ConductancePriorityQueue is a MaxHeap => binary tree
    if (SuperPQ::size() > 2) {
      f = std::max(f, SuperPQ::heap[2].value);
    }
    unlock(syncronized);
    return f;
  }

  double_t topConductance(bool syncronized = false) const {
    ConductanceFraction f = topFraction(synchronized);
    return f.value();
  }
  double_t secondTopCondunctance(bool syncronized = false) const {
    ConductanceFraction f = secondTopFraction(synchronized);
    return f.value();
  }

  // ! Get the conductance fraction of a partition
  // (better use PartitionedHypergraph getPartCutWeight and getPartVolume)
  ConductanceFraction getFraction(const PartitionID p, bool syncronized = false) const {
    lock(synchronized);
    ConductanceFraction f = SuperPQ::getKey(p);
    unlock(synchronized);
    return f;
  }
  // ! Get the conductance of a partition
  // (better use PartitionedHypergraph getPartCutWeight and getPartVolume)
  double_t getConductance(const PartitionID p, bool synchronized = false) const {
    lock(synchronized);
    ConductanceFraction f = SuperPQ::getKey(p);
    unlock(synchronized);
    return f.value();
  }

  // ################# MUTEX OPERATIONS FOR CHANGING PQ #################

  void lock(bool synchronized) {
    if synchronized _pq_lock.lock();
  }

  void unlock(bool synchronized) {
    if synchronized _pq_lock.unlock();
  }

  // ################# MEMBER VARIABLES #################
  SpinLock _pq_lock;
  HyperedgeWeight _total_volume;
  PartitionID _size;
  vec<bool> _complement_val_bits;
  bool _initialized;
};


}  // namespace ds
}  // namespace mt_kahypar