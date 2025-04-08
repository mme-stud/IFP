
#pragma once

#include <atomic>
#include <mutex>

#include <tbb/parallel_invoke.h>

#include "kahypar-resources/meta/mandatory.h"

#include "mt-kahypar/datastructures/partitioned_hypergraph.h"
#include "mt-kahypar/datastructures/priority_queue.h"
#include "mt-kahypar/datastructures/nonnegative_fraction.h"

#include "mt-kahypar/datastructures/delta_val.h"

namespace mt_kahypar {
namespace ds {
  
/* Moved to hypergraph_common.h

using ConductanceFraction = NonnegativeFraction<HypergraphVolume>;

struct ConductanceInfo {
  ConductanceFraction fraction;
  PartitionID partID;
};
*/


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
  using DeltaV = DeltaValue<HypergraphVolume>;
public:
  ConductancePriorityQueue() :
    SuperPQ(0),
    _total_volume(-1),
    _size(0),
    _complement_val_bits(),
    _delta_part_volumes(),
    _delta_cut_weights(),
    _initialized(false)
    { }
  
  // ! Initializes the priority queue with the partitions of the hypergraph
  // ! Could be called concurrently !!!
  void initialize(const PartitionedHypergraph& hg, bool synchronized = false) {
    /// [debug] std::cerr << "ConductancePriorityQueue::initialize(hg, " << V(synchronized) << ")" << std::endl;
    // ASSERT(!_initialized); could be called concurrently -> no assertion before lock
    lock(synchronized);
    if (_initialized) {
      unlock(synchronized);
      return;
    }
    _total_volume = getHGTotalVolume(hg);
    _size = hg.k();
    _complement_val_bits.resize(_size);
    _delta_part_volumes.resize(_size, 0);
    _delta_cut_weights.resize(_size, 0);
    SuperPQ::clear();
    SuperPQ::resize(_size);
    SuperPQ::heap.resize(_size);
    tbb::parallel_for(PartitionID(0), _size, [&](const PartitionID& p) {
      HypergraphVolume cut_weight = getHGPartCutWeight(hg, p);
      HypergraphVolume part_volume = getHGPartVolume(hg, p);
      _complement_val_bits[p] = (part_volume > _total_volume - part_volume);
      ConductanceFraction f(cut_weight, std::min(part_volume, _total_volume - part_volume));
      SuperPQ::heap[p].id = p; 
      SuperPQ::heap[p].key = f;
      SuperPQ::positions[p] = p;
    });
    buildHeap();
    _initialized = true;
    unlock(synchronized);
  }

  // ! Returns true if the priority queue is initialized
  bool initialized() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::initialized()" << std::endl;
    return _initialized;
  }

  // ! Reset the priority queue to the uninitialized state
  void reset(bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::reset()" << std::endl;
    lock(synchronized);
    SuperPQ::clear();
    _total_volume = -1;
    _size = 0;
    _complement_val_bits.clear();
    _delta_part_volumes.clear();
    _delta_cut_weights.clear();
    _initialized = false;
    unlock(synchronized);
  }

  // ! Returns an approximate memory consumption of the conductance priority queue in bytes
  size_t memoryConsumption() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::memoryConsumption()" << std::endl;
    return SuperPQ::memoryConsumption() + _complement_val_bits.size() * sizeof(bool)
    + (_delta_cut_weights.size() + _delta_part_volumes.size()) * (sizeof(bool) + sizeof(HypergraphVolume));
  }

  // ! Updates the priority queue after global changes in partition
  void globalUpdate(const PartitionedHypergraph& hg, bool synchronized = false) {
    /// [debug] std::cerr << "ConductancePriorityQueue::globalUpdate(hg, " << V(synchronized) << ")" << std::endl;
    lock(synchronized);
    ASSERT(_initialized && _size == hg.k() && _size == static_cast<PartitionID>(size()));
    if (_uses_original_stats) {
      ASSERT(_total_volume == hg.originalTotalVolume(), "Total volume in ConductancePriorityQueue is" << _total_volume << ", but should be" << hg.originalTotalVolume());
    }
    _total_volume = getHGTotalVolume(hg);
    tbb::parallel_for(PartitionID(0), _size, [&](const PartitionID& p) {
      HypergraphVolume cut_weight = getHGPartCutWeight(hg, p);
      HypergraphVolume part_volume = getHGPartVolume(hg, p);
      ASSERT(part_volume <= _total_volume, "Partition volume " << part_volume << " is greater than total volume " << _total_volume);
      ASSERT(cut_weight <= part_volume, "Cut weight " << cut_weight << " is greater than partition volume " << part_volume);
      ASSERT(cut_weight + part_volume <= _total_volume, "Cut weight " << cut_weight << " plus partition volume " << part_volume << " is greater than total volume " << _total_volume);
      _complement_val_bits[p] = (part_volume > _total_volume - part_volume);
      ConductanceFraction f(cut_weight, std::min(part_volume, _total_volume - part_volume));
      SuperPQ::heap[SuperPQ::positions[p]].key = f;
      ASSERT(_delta_cut_weights[p] == 0 && _delta_part_volumes[p] == 0, "Deltas should be empty: " << V(_delta_cut_weights[p]) << ", " V(_delta_part_volumes[p]));
    });
    buildHeap();
    unlock(synchronized);
  }

  // ! Checks if the priority queue is correct with respect to the hypergraph: const version
  bool check(const PartitionedHypergraph& hg) const {
    /// [debug] std::cerr << "ConductancePriorityQueue::check(hg)" << std::endl;
    ASSERT(_initialized && _size == hg.k());
    bool correct = true;
    if (_total_volume != getHGTotalVolume(hg)) {
      correct = false;
      LOG << "Total volume in ConductancePriorityQueue is" << _total_volume << ", but should be" << getHGTotalVolume(hg);
    }
    for (PartitionID p = 0; p < _size; ++p) {
      ConductanceFraction f = SuperPQ::getKey(p);
      HypergraphVolume cut_weight = f.getNumerator();
      HypergraphVolume part_volume = f.getDenominator();
      ASSERT(part_volume <= _total_volume);
      ASSERT(cut_weight <= part_volume);
      if (_complement_val_bits[p]) {
        part_volume = _total_volume - part_volume;
      }
      ASSERT(cut_weight <= part_volume);
      if (part_volume != getHGPartVolume(hg, p)) {
        correct = false;
        LOG << "Volume of partition in ConductancePriorityQueue" << V(p) << "is" << V(part_volume) << ", but should be" << getHGPartVolume(hg, p);
      }
      if (cut_weight != getHGPartCutWeight(hg, p)) {
        correct = false;
        LOG << "Cut weight of partition in ConductancePriorityQueue" << V(p) << "is" << V(cut_weight) << ", but should be" << getHGPartCutWeight(hg, p);
      }
      // Deltas should be empty after changeNodePart-stage is empty. And check is called only then
      if (_delta_part_volumes[p] != 0 || _delta_cut_weights[p] != 0) {
        correct = false;
        LOG << "Deltas should be always empty, when changeNodePart-stage is finished: " << V(_delta_part_volumes[p]) << ", " << V(_delta_cut_weights[p]);
      }
    }
    correct = correct && SuperPQ::isHeap() && SuperPQ::positionsMatch();
    return correct;
  }

  // ! Checks if the priority queue is correct with respect to the hypergraph: synchronizable version
  bool checkSync(const PartitionedHypergraph& hg, bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::checkSync(hg, " << V(synchronized) << ")" << std::endl;
    lock(synchronized);
    bool correct = check(hg);
    unlock(synchronized);
    return correct;
  }

  // ################# Priority Queue Operations #################

  // used to upderstand, if buildHeap is nessesary
  // otherwise use SuperPQ::isHeap() - it prints out to LOG
  bool isHeap() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::isHeap()" << std::endl;
    for (PartitionID i = 1; i < size(); ++i) {
      if (heap[parent(i)].key < heap[i].key) {
        // LOG << "heap property violation" << V(i) << V(parent(i))  << V(heap[i].key) << V(heap[parent(i)].key);
        return false;
      }
    }
    return true;
  }

  // ! Adjusts the cut weight and the volume of a partition by values
  // ! changes pq => uses a lock  
  // ! Needs exact new values => discouraged from usage due to patential data race 
  void adjustKey(const PartitionID& p, const HypergraphVolume& cut_weight, const HypergraphVolume& part_volume, bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::adjustKey(" << V(p) << ", " << V(cut_weight) << ", " << V(part_volume) << ", " << V(synchronized) << ")" << std::endl;
    ASSERT(_initialized);
    ASSERT(static_cast<size_t>(_size) == _complement_val_bits.size());
    if (part_volume < cut_weight || _total_volume < part_volume || _total_volume < part_volume + cut_weight) {
      // this update is incorrect (potentially due to concurrency) => will be redone later
      // [used by changeNodePart, where only for the last thread changing partition p 
      //                                  the right stats are guaranteed]
      // LOG << "ConductancePriorityQueue::adjustKey(" << V(p) << ", " << V(cut_weight) << ", " << V(part_volume) << ") is skipped due to incorrect stats: " << " " << V(_total_volume) << ". " << "[shouldn't be a problem]";
      return;
    }
    lock(synchronized);
    // LOG << "ConductancePriorityQueue::adjustKey(" << V(p) << ", " << V(cut_weight) << ", " << V(part_volume) << ") is started";
    _complement_val_bits[p] = (part_volume > _total_volume - part_volume);
    ConductanceFraction f(cut_weight, std::min(part_volume, _total_volume - part_volume));
    SuperPQ::adjustKey(p, f);
    SuperPQ::heap[SuperPQ::positions[p]].key = f; // needed, as the fraction are equal if their reduced forms are equal
    // so SuperPQ::adjustKey(p, f) would not change the key to the new one
    // but in ConductancePriorityQueue we need to set the key to the exact numerator and denominator
    // as they have meaning in the context of the hypergraph
    // LOG << "ConductancePriorityQueue::adjustKey(" << V(p) << ", " << V(cut_weight) << ", " << V(part_volume) << ") is finished";
    unlock(synchronized);
  }
  
  // ! Adjusts the cut weight and the volume of a partition by deltas of theis values
  // ! Used only by changeNodePart, where only for the last thread changing partition p 
  // ! the right stats are guaranteed
  // changes pq => uses a lock
  // Skips clearly not finished updates (_total_value < _part_volume etc)
  // Maintains deltas of skipped updates
  void adjustKeyByDeltas(const PartitionID& p, 
                         const DeltaV& d_cut_weight, 
                         const DeltaV& d_part_volume, 
                         bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::adjustKeyByDeltas(" << V(p) << ", " << V(d_cut_weight) << ", " << V(d_part_volume) << ", " << V(synchronized) << ")" << std::endl;
    ASSERT(_initialized);
    ASSERT(static_cast<size_t>(_size) == _complement_val_bits.size());
    lock(synchronized);
    // Update daltas
    _delta_part_volumes[p] += d_part_volume;
    _delta_cut_weights[p] += d_cut_weight;

    // Get current values
    ConductanceFraction f = SuperPQ::getKey(p);
    HypergraphVolume cut_weight = f.getNumerator();
    HypergraphVolume part_volume = f.getDenominator();
    if (_complement_val_bits[p]) {
      ASSERT(_total_volume >= part_volume);
      part_volume = _total_volume - part_volume;
    }
    // Check if update could be skipped as clearly not finished
    DeltaV new_cut_weight = _delta_cut_weights[p] + cut_weight;
    if (new_cut_weight.isNegative() || new_cut_weight.abs() > _total_volume) {
      // skip this update, as it is not finished
      // LOG << "ConductancePriorityQueue::adjustKeyByDeltas(" << V(p) << ", " << V(d_cut_weight) << ", " << V(d_part_volume) << ") is skipped due to incorrect stats: " << " " << V(_total_volume) << ", " << V(new_cut_weight) << ", " << V(new_part_volume) << ". [shouldn't be a problem]";
      unlock(synchronized);
      return;
    }
    DeltaV new_part_volume = _delta_part_volumes[p] + part_volume;
    if (new_part_volume.isNegative() || new_part_volume.abs() > _total_volume || (new_cut_weight + new_part_volume).abs() > _total_volume) {
      // skip this update, as it is not finished
      // LOG << "ConductancePriorityQueue::adjustKeyByDeltas(" << V(p) << ", " << V(d_cut_weight) << ", " << V(d_part_volume) << ") is skipped due to incorrect stats: " << " " << V(_total_volume) << ", " << V(new_cut_weight) << ", " << V(new_part_volume) << ". [shouldn't be a problem]";
      unlock(synchronized);
      return;
    }
   
    // Didn't skip the update
    // LOG << "ConductancePriorityQueue::adjustKeyByDeltas(" << V(p) << ", " << V(new_cut_weight) << ", " << V(new_part_volume) << ") is started";
    
    // Clear deltas
    _delta_cut_weights[p] = 0;  _delta_part_volumes[p] = 0;
    // Adjust key
    _complement_val_bits[p] = (new_part_volume.abs() > _total_volume - new_part_volume.abs());
    ConductanceFraction newF(new_cut_weight.abs(), std::min(new_part_volume.abs(), _total_volume - new_part_volume.abs()));
    SuperPQ::adjustKey(p, newF);
    SuperPQ::heap[SuperPQ::positions[p]].key = newF; // needed, as the fraction are equal if their reduced forms are equal
    // so SuperPQ::adjustKey(p, f) would not change the key to the new one
    // but in ConductancePriorityQueue we need to set the key to the exact numerator and denominator
    // as they have meaning in the context of the hypergraph

    // LOG << "ConductancePriorityQueue::adjustKeyByDeltas(" << V(p) << ", " << V(new_cut_weight) << ", " << V(new_part_volume) << ") is finished";
    unlock(synchronized);
  }

  // ! Updates PQ after total volume of the hypergraph has changed
  // ! changes pq => uses a lock
  void updateTotalVolume(const HypergraphVolume& new_total_volume, bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::updateTotalVolume(" << V(new_total_volume) << ", " << V(synchronized) << ")" << std::endl;		
    ASSERT(_initialized && static_cast<size_t>(_size) == _complement_val_bits.size());
    lock(synchronized);
    for (PartitionID p = 0; p < _size; ++p) {
      ConductanceFraction& f = SuperPQ::getKey(p);
      HypergraphVolume part_volume = f.getDenominator();
      ASSERT(part_volume <= _total_volume && part_volume <= new_total_volume);
      if (_complement_val_bits[p]) {
        part_volume = _total_volume - part_volume;
      }
      _complement_val_bits[p] = (part_volume > new_total_volume - part_volume);
      f.setDenominator(std::min(part_volume, new_total_volume - part_volume));
      ASSERT(_delta_cut_weights[p] == 0 && _delta_part_volumes[p] == 0, "Deltas should be empty: " << V(_delta_cut_weights[p]) << ", " V(_delta_part_volumes[p]));
    }
    buildHeap();
    _total_volume = new_total_volume;
    unlock(synchronized);
  }

  // ! Get the partition with the highest conductance: const version
  ConductanceInfo top() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::top()" << std::endl;
    auto first = SuperPQ::heap[0];
    return ConductanceInfo(first.key, first.id);
  }

  // ! Get the partition with the highest conductance: synchronizable version
  ConductanceInfo topSync(bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::topSync(" << V(synchronized) << ")" << std::endl;
    lock(synchronized);
    ConductanceInfo first = top();
    unlock(synchronized);
    return first;
  }

  // ! Get the partition with the second highest conductance: const version
  ConductanceInfo secondTop() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::secondTop()" << std::endl;
    ASSERT(SuperPQ::size() > 1);
    auto second = SuperPQ::heap[1];
    // ConductancePriorityQueue is a MaxHeap => binary tree
    if (size() > 2 && SuperPQ::heap[1].key < SuperPQ::heap[2].key) {
      second = SuperPQ::heap[2];
    }
    return ConductanceInfo(second.key, second.id);
  }
  
  // ! Get the partition with the second highest conductance: synchronizable version
  ConductanceInfo secondTopSync(bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::secondTopSync(" << V(synchronized) << ")" << std::endl;
    lock(synchronized);
    ConductanceInfo second = secondTop();
    unlock(synchronized);
    return second;
  }

  // ! Get ConductanceInfo of the top three partitions (unsorted, [0] is max): constant version
  // ! (Works only for a binary heap)
  vec<ConductanceInfo> topThree() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::topThree()" << std::endl;
    vec<ConductanceInfo> top_three(3, ConductanceInfo());
    for (PartitionID i = 0; i < 3 && i < _size; ++i) {
      auto elem = SuperPQ::heap[i]; // HeapElement { KeyT, IdT }
      top_three[i] = ConductanceInfo(elem.key, elem.id);
    }
    return top_three;
  }

  // ! Get ConductanceInfo of the top three partitions (unsorted, [0] is max): synchronizable version
  // ! (Works only for a binary heap)
  vec<ConductanceInfo> topThreeSync(bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::topThreeSync(" << V(synchronized) << ")" << std::endl;
    lock(synchronized);
    vec<ConductanceInfo> top_three = topThree();
    unlock(synchronized);
    return top_three;
  }
  
  bool empty() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::empty()" << std::endl;
    bool e = SuperPQ::empty();
    return e;
  }

  PartitionID size() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::size()" << std::endl;
    ASSERT(static_cast<PosT>(_size) == SuperPQ::size());
    return _size;
  }

  // ################## USAGE OF ORIGINAL PHG STATS #################

  // ! Uses the original stats of the hypergraph
  bool usesOriginalStats() const {
    /// [debug] std::cerr << "ConductancePriorityQueue::usesOriginalStats()" << std::endl;
    return _uses_original_stats;
  }
  
  // ! Makes ConductancePriorityQueue use current _total_volume and _part_volumes
  // ! instead of the original (inherited from old versions of hg) ones
  // ! To be used before initialization
  void disableUsageOfOriginalHGStats() {
    /// [debug] std::cerr << "ConductancePriorityQueue::disableUsageOfOriginalHGStats()" << std::endl;
    ASSERT(!_initialized, "ConductancePriorityQueue is already initialized");
    _uses_original_stats = false;
  }

  // ! Makes ConductancePriorityQueue use original _total_volume and _part_volumes
  // ! instead of the current ones
  // ! To be used before initialization
  void enableUsageOfOriginalHGStats() {
    /// [debug] std::cerr << "ConductancePriorityQueue::enableUsageOfOriginalHGStats()" << std::endl;
    ASSERT(!_initialized, "ConductancePriorityQueue is already initialized");
    _uses_original_stats = true;
  }

  // ################# MUTEX OPERATIONS FOR CHANGING PQ #################

  // ! Exclude simultaneous synchronized access to the priority queue
  // ! To be used when HG is changed and the priority queue is changed in parallel
  void lock(bool synchronized = true) {
    if (synchronized) _pq_lock.lock();
      /// [debug] std::cerr << "ConductancePriorityQueue::lock(" << V(synchronized) << ")" << std::endl;
  }

  void unlock(bool synchronized = true) {
    if (synchronized) { 
      /// [debug] std::cerr << "ConductancePriorityQueue::unlock(" << V(synchronized) << ")" << std::endl;
      ASSERT(!_pq_lock.tryLock(), "ConductancePriorityQueue::unlock() called without lock");
      _pq_lock.unlock();
    }
  }

private:
  // ! Builds the heap in O(_size) time
  // ! no built in lock
  void buildHeap() {
    /// [debug] std::cerr << "ConductancePriorityQueue::buildHeap()" << std::endl;
    ASSERT(static_cast<PosT>(_size) == SuperPQ::size());
    if (!isHeap()) {
      for (PartitionID p = _size - 1; p >= 0; --p) {
        SuperPQ::siftDown(p);
      }
    }
    ASSERT(SuperPQ::isHeap() && SuperPQ::positionsMatch());
  }

  // ################### COMMUNICATION WITH THE HG ######################
  // ! Get needed kind of total volume
  // ! (original or current)
  HypergraphVolume getHGTotalVolume(const PartitionedHypergraph& hg) const {
    /// [debug] std::cerr << "ConductancePriorityQueue::getHGTotalVolume(hg): " << hg.originalTotalVolume() << ", " << hg.totalVolume() << std::endl;
    if (_uses_original_stats) {
      return hg.originalTotalVolume();
    } else {
      return hg.totalVolume();
    }
  }
  
  // ! Get needed kind of part volume
  // ! (original or current)
  HypergraphVolume getHGPartVolume(const PartitionedHypergraph& hg, const PartitionID p) const {
    /// [debug] std::cerr << "ConductancePriorityQueue::getHGPartVolume(hg, " << p << "): "  << hg.partOriginalVolume(p) << ", " << hg.partVolume(p) << std::endl;
    if (_uses_original_stats) {
      return hg.partOriginalVolume(p);
    } else {
      return hg.partVolume(p);
    }
  }

  // ! Get needed kind of cut weight
  // ! (only one kind of cut weight for now)
  HypergraphVolume getHGPartCutWeight(const PartitionedHypergraph& hg, const PartitionID p) const {
    /// [debug] std::cerr << "ConductancePriorityQueue::getHGPartCutWeight(hg, " << p << "): " << hg.partCutWeight(p) << std::endl;
    return hg.partCutWeight(p);
  }

  // #################### POTENTIALLY USELESS PART ######################
  
  // ! insert a partition with its cut weight and volume
  // ! changes pq => uses a lock
  // ! not sure if this is useful
  void insert(const PartitionID& p, const HypergraphVolume& cut_weight, const HypergraphVolume& volume, bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::insert(" << p << ", " << cut_weight << ", " << volume << ", " << synchronized << ")" << std::endl;
    ConductanceFraction f(cut_weight, std::min(volume, _total_volume - volume));
    lock(synchronized);
    ASSERT(_total_volume >= volume);
    _complement_val_bits[p] = (volume > _total_volume - volume);
    SuperPQ::insert(p, f);
    _size++;
    unlock(synchronized);
  }

  // ! changes pq => uses a lock
  // ! not sure if this is useful
  void remove(const PartitionID& p, bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::remove(" << p << ", " << synchronized << ")" << std::endl;
    lock(synchronized);
    SuperPQ::remove(p);
    _size--;
    unlock(synchronized);
  }


  // ! not sure if this is useful
  void deleteTop(bool synchronized = true) {
    /// [debug] std::cerr << "ConductancePriorityQueue::deleteTop(" << synchronized << ")" << std::endl;
    lock(synchronized);
    SuperPQ::deleteTop();
    _size--;
    unlock(synchronized);
  }

  // ! Get the conductance fraction of the partition with the highest conductance
  ConductanceFraction topFractionSync(bool synchronized = false) {
    /// [debug] std::cerr << "ConductancePriorityQueue::topFractionSync(" << synchronized << ")" << std::endl;
    lock(synchronized);
    ConductanceFraction f = SuperPQ::topKey();
    unlock(synchronized);
    return f;
  }
  // ! Get the conductance fraction of the partition with the second highest conductance
  ConductanceFraction secondTopFractionSync(bool synchronized = false) {
    /// [debug] std::cerr << "ConductancePriorityQueue::secondTopFractionSync(" << synchronized << ")" << std::endl;
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
    /// [debug] std::cerr << "ConductancePriorityQueue::topConductanceSync(" << synchronized << ")" << std::endl;
    ConductanceFraction f = topFractionSync(synchronized);
    return f.value();
  }
  double_t secondTopCondunctanceSync(bool synchronized = false) {
    // [debug] std::cerr << "ConductancePriorityQueue::secondTopConductanceSync(" << synchronized << ")" << std::endl;
    ConductanceFraction f = secondTopFractionSync(synchronized);
    return f.value();
  }

  // ! Get the conductance fraction of a partition
  // (better use PartitionedHypergraph getPartCutWeight and getPartVolume)
  ConductanceFraction getFractionSync(const PartitionID p, bool synchronized = false) {
    /// [debug] std::cerr << "ConductancePriorityQueue::getFractionSync(" << p << ", " << synchronized << ")" << std::endl;
    lock(synchronized);
    ConductanceFraction f = SuperPQ::getKey(p);
    unlock(synchronized);
    return f;
  }
  // ! Get the conductance of a partition
  // (better use PartitionedHypergraph getPartCutWeight and getPartVolume)
  double_t getConductanceSync(const PartitionID p, bool synchronized = false) {
    /// [debug] std::cerr << "ConductancePriorityQueue::getConductanceSync(" << p << ", " << synchronized << ")" << std::endl;
    lock(synchronized);
    ConductanceFraction f = SuperPQ::getKey(p);
    unlock(synchronized);
    return f.value();
  }

  // ################# MEMBER VARIABLES #################
  SpinLock _pq_lock;
  HypergraphVolume _total_volume;
  PartitionID _size;
  vec<bool> _complement_val_bits;
  vec<DeltaV> _delta_part_volumes;
  vec<DeltaV> _delta_cut_weights;
  bool _uses_original_stats = true;
  bool _initialized;
};


}  // namespace ds
}  // namespace mt_kahypar