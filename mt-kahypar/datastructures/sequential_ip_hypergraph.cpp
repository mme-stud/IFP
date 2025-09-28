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

#include "mt-kahypar/datastructures/sequential_ip_hypergraph.h"
// #include "mt-kahypar/datastructures/concurrent_bucket_map.h"
#include "mt-kahypar/datastructures/sequential_bucket_map.h"
#include "mt-kahypar/utils/timer.h"
#include "mt-kahypar/utils/memory_tree.h"


namespace mt_kahypar::ds {

  /*!
  * This struct is used during multilevel coarsening to efficiently
  * detect parallel hyperedges.
  */
  struct ContractedHyperedgeInformation {
    HyperedgeID he = kInvalidHyperedge;
    size_t hash = kEdgeHashSeed;
    size_t size = std::numeric_limits<size_t>::max();
    bool valid = false;
  };

   /*!
   * Contracts the community structure given in `_part_ids`. All vertices with
   * the same label are collapsed into the same vertex. The resulting single-pin
   * (if single-pin nets removal isn't disabled) and parallel hyperedges are
   * removed from the contracted graph. The function maps the passed community
   * structure to the new community labels in the coarse hypergraph.
   *
   * \param global_communities Community structure that should be mapped to the 
   *                           new community labels in the coarse hypergraph.
   */
  void SequentialIPHypergraph::contract(std::vector<PartitionID>& global_communities, bool deterministic) {

    if ( !_tmp_contraction_buffer ) {
      allocateTmpContractionBuffer();
    }

    // Auxiliary buffers - reused during multilevel hierarchy to prevent expensive allocations
    Array<size_t>& mapping = _tmp_contraction_buffer->mapping;
    Array<Hypernode>& tmp_hypernodes = _tmp_contraction_buffer->tmp_hypernodes;
    IncidentNets& tmp_incident_nets = _tmp_contraction_buffer->tmp_incident_nets;
    Array<size_t>& tmp_num_incident_nets =
            _tmp_contraction_buffer->tmp_num_incident_nets;
    Array<HypernodeWeight>& hn_weights =
            _tmp_contraction_buffer->hn_weights;
    Array<Hyperedge>& tmp_hyperedges = _tmp_contraction_buffer->tmp_hyperedges;
    IncidenceArray& tmp_incidence_array = _tmp_contraction_buffer->tmp_incidence_array;
    Array<size_t>& he_sizes = _tmp_contraction_buffer->he_sizes;
    Array<size_t>& valid_hyperedges = _tmp_contraction_buffer->valid_hyperedges;

    ASSERT(static_cast<size_t>(_num_hypernodes) <= mapping.size());
    ASSERT(static_cast<size_t>(_num_hypernodes) <= tmp_hypernodes.size());
    ASSERT(static_cast<size_t>(_total_degree) <= tmp_incident_nets.size());
    ASSERT(static_cast<size_t>(_num_hypernodes) <= tmp_num_incident_nets.size());
    ASSERT(static_cast<size_t>(_num_hypernodes) <= hn_weights.size());
    ASSERT(static_cast<size_t>(_num_hyperedges) <= tmp_hyperedges.size());
    ASSERT(static_cast<size_t>(_num_pins) <= tmp_incidence_array.size());
    ASSERT(static_cast<size_t>(_num_hyperedges) <= he_sizes.size());
    ASSERT(static_cast<size_t>(_num_hyperedges) <= valid_hyperedges.size());


    // #################### STAGE 1 ####################
    // Compute vertex ids of coarse hypergraph with a parallel prefix sum
    mapping.assign(_num_hypernodes, 0);
    std::vector<PartitionID> communities = _part_ids;

    doForAllNodes([&](const HypernodeID& hn) {
      ASSERT(static_cast<size_t>(communities[hn]) < mapping.size());
      mapping[communities[hn]] = UL(1);
    });


    // Prefix sum determines vertex ids in coarse hypergraph
    std::vector<size_t> mapping_prefix_sum(_num_hypernodes + 1);
    mapping_prefix_sum[0] = 0;
    for ( size_t i = 0; i < UL(_num_hypernodes); ++i ) {
      mapping_prefix_sum[i + 1] = mapping[i] + mapping_prefix_sum[i];
    }
    HypernodeID num_hypernodes = mapping_prefix_sum[_num_hypernodes];

    // Remap community ids
    for ( HypernodeID hn = 0; hn < _num_hypernodes; ++hn ) {
      if ( nodeIsEnabled(hn) ) {
        communities[hn] = mapping_prefix_sum[communities[hn]];
      } else {
        communities[hn] = kInvalidHypernode;
      }

      // Reset tmp contraction buffer
      if ( hn < num_hypernodes ) {
        hn_weights[hn] = 0;
        tmp_hypernodes[hn] = Hypernode(true);
        tmp_num_incident_nets[hn] = 0;
      }
    }

    // Mapping from a vertex id of the current hypergraph to its
    // id in the coarse hypergraph
    auto map_to_coarse_hypergraph = [&](const HypernodeID hn) {
      ASSERT(hn < communities.size());
      return communities[hn];
    };

    doForAllNodes([&](const HypernodeID& hn) {
      const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
      ASSERT(coarse_hn < num_hypernodes, V(coarse_hn) << V(num_hypernodes));
      // Weight vector is atomic => thread-safe
      hn_weights[coarse_hn] += nodeWeight(hn);
      // Aggregate upper bound for number of incident nets of the contracted vertex
      tmp_num_incident_nets[coarse_hn] += nodeDegree(hn);
    });

    // #################### STAGE 2 ####################
    // In this step hyperedges and incident nets of vertices are contracted inside the temporary
    // buffers. The vertex ids of pins are already remapped to the vertex ids in the coarse
    // graph and duplicates are removed. Also nets that become single-pin hyperedges are marked
    // as invalid. All incident nets of vertices that are collapsed into one vertex in the coarse
    // graph are also aggregate in a consecutive memory range and duplicates are removed. Note
    // that parallel and single-pin hyperedges are not removed from the incident nets (will be done
    // in a postprocessing step).
    // !!! (New) If single-pin nets removal is disabled, we do not invalidate single-pin hyperedges
    
    auto cs2 = [](const HypernodeID x) { return x * x; };
    // ConcurrentBucketMap<ContractedHyperedgeInformation> hyperedge_hash_map;
    SequentialBucketMap<ContractedHyperedgeInformation> hyperedge_hash_map;
    hyperedge_hash_map.reserve_for_estimated_number_of_insertions(_num_hyperedges);
    
    // Contract Hyperedges
    {
      for (HyperedgeID he = HyperedgeID(0); he < _num_hyperedges; ++he) {
        if ( edgeIsEnabled(he) ) {
          // Copy hyperedge and pins to temporary buffer
          const Hyperedge& e = _hyperedges[he];
          ASSERT(static_cast<size_t>(he) < tmp_hyperedges.size());
          ASSERT(e.firstInvalidEntry() <= tmp_incidence_array.size());
          tmp_hyperedges[he] = e;
          valid_hyperedges[he] = 1;

          // Map pins to vertex ids in coarse graph
          const size_t incidence_array_start = tmp_hyperedges[he].firstEntry();
          const size_t incidence_array_end = tmp_hyperedges[he].firstInvalidEntry();
          for ( size_t pos = incidence_array_start; pos < incidence_array_end; ++pos ) {
            const HypernodeID pin = _incidence_array[pos];
            ASSERT(pos < tmp_incidence_array.size());
            tmp_incidence_array[pos] = map_to_coarse_hypergraph(pin);
          }

          // Remove duplicates and disabled vertices
          auto first_entry_it = tmp_incidence_array.begin() + incidence_array_start;
          std::sort(first_entry_it, tmp_incidence_array.begin() + incidence_array_end);
          auto first_invalid_entry_it = std::unique(first_entry_it, tmp_incidence_array.begin() + incidence_array_end);
          while ( first_entry_it != first_invalid_entry_it && *(first_invalid_entry_it - 1) == kInvalidHypernode ) {
            --first_invalid_entry_it;
          }

          // Update size of hyperedge in temporary hyperedge buffer
          const size_t contracted_size = std::distance(
                  tmp_incidence_array.begin() + incidence_array_start, first_invalid_entry_it);
          tmp_hyperedges[he].setSize(contracted_size);


          if ( contracted_size > 1 || _disable_single_pin_nets_removal) {
            // Compute hash of contracted hyperedge
            size_t footprint = kEdgeHashSeed;
            for ( size_t pos = incidence_array_start; pos < incidence_array_start + contracted_size; ++pos ) {
              footprint += cs2(tmp_incidence_array[pos]);
            }
            hyperedge_hash_map.insert(footprint,
                                      ContractedHyperedgeInformation{ he, footprint, contracted_size, true });
          } else {
            // Hyperedge becomes a single-pin hyperedge
            // single-pin nets removal is not disabled
            valid_hyperedges[he] = 0;
            tmp_hyperedges[he].disable();
          }
        } else {
          valid_hyperedges[he] = 0;
        }
      }
    }

    // Contract Incident Nets
    // Compute start position the incident nets of a coarse vertex in the
    // temporary incident nets array with a parallel prefix sum
    {
      std::vector<size_t> tmp_incident_nets_pos(num_hypernodes, 0);
      std::vector<size_t> tmp_incident_nets_prefix_sum(num_hypernodes + 1, 0);
      tmp_incident_nets_prefix_sum[0] = 0;
      for ( size_t i = 0; i < UL(num_hypernodes); ++i ) {
        tmp_incident_nets_prefix_sum[i + 1] = tmp_num_incident_nets[i] + tmp_incident_nets_prefix_sum[i];
      }

      // Write the incident nets of each contracted vertex to the temporary incident net array
      doForAllNodes([&](const HypernodeID& hn) {
        const HypernodeID coarse_hn = map_to_coarse_hypergraph(hn);
        const HyperedgeID node_degree = nodeDegree(hn);
        size_t incident_nets_pos = tmp_incident_nets_prefix_sum[coarse_hn] +
                                   tmp_incident_nets_pos[coarse_hn];
        tmp_incident_nets_pos[coarse_hn] += node_degree;
        ASSERT(incident_nets_pos + node_degree <= tmp_incident_nets_prefix_sum[coarse_hn + 1]);
        memcpy(tmp_incident_nets.data() + incident_nets_pos,
               _incident_nets.data() + _hypernodes[hn].firstEntry(),
               sizeof(HyperedgeID) * node_degree);
      });

      
      // Setup temporary hypernodes
      std::vector<HypernodeID> high_degree_vertices;
      for (HypernodeID coarse_hn = 0; coarse_hn < num_hypernodes; ++coarse_hn) {
        // Remove duplicates
        const size_t incident_nets_start = tmp_incident_nets_prefix_sum[coarse_hn];
        const size_t incident_nets_end = tmp_incident_nets_prefix_sum[coarse_hn + 1];
        const size_t tmp_degree = incident_nets_end - incident_nets_start;
        if ( tmp_degree <= HIGH_DEGREE_CONTRACTION_THRESHOLD ) {
          std::sort(tmp_incident_nets.begin() + incident_nets_start,
                    tmp_incident_nets.begin() + incident_nets_end);
          auto first_invalid_entry_it = std::unique(tmp_incident_nets.begin() + incident_nets_start,
                                                    tmp_incident_nets.begin() + incident_nets_end);
          
          // Setup pointers to temporary incident nets
          const size_t contracted_size = std::distance(tmp_incident_nets.begin() + incident_nets_start,
                                                       first_invalid_entry_it);
          tmp_hypernodes[coarse_hn].setSize(contracted_size);
        } else {
          high_degree_vertices.push_back(coarse_hn);
        }
        tmp_hypernodes[coarse_hn].setWeight(hn_weights[coarse_hn]);
        tmp_hypernodes[coarse_hn].setFirstEntry(incident_nets_start);
      }
      
      if ( !high_degree_vertices.empty() ) {
        // High degree vertices are treated special, because sorting and afterwards
        // removing duplicates can become a major sequential bottleneck. Therefore,
        // we distribute the incident nets of a high degree vertex into our sequential
        // bucket map. As a result all equal incident nets reside in the same bucket
        // afterwards. In a second step, we process each bucket and apply
        // for each bucket the duplicate removal procedure from above.
        SequentialBucketMap<HyperedgeID> duplicate_incident_nets_map;
        // ConcurrentBucketMap<HyperedgeID> duplicate_incident_nets_map;
        for ( const HypernodeID& coarse_hn : high_degree_vertices ) {
          const size_t incident_nets_start = tmp_incident_nets_prefix_sum[coarse_hn];
          const size_t incident_nets_end = tmp_incident_nets_prefix_sum[coarse_hn + 1];
          const size_t tmp_degree = incident_nets_end - incident_nets_start;
          
          // Insert incident nets into concurrent bucket map
          duplicate_incident_nets_map.reserve_for_estimated_number_of_insertions(tmp_degree);
          for(size_t pos = incident_nets_start; pos < incident_nets_end; ++pos) {
            HyperedgeID he = tmp_incident_nets[pos];
            duplicate_incident_nets_map.insert(he, std::move(he));
          }

          // Process each bucket and remove duplicates
          size_t incident_nets_pos = incident_nets_start;
          for (size_t bucket = 0; bucket < duplicate_incident_nets_map.numBuckets(); ++bucket) {
            auto& incident_net_bucket = duplicate_incident_nets_map.getBucket(bucket);
            std::sort(incident_net_bucket.begin(), incident_net_bucket.end());
            auto first_invalid_entry_it = std::unique(incident_net_bucket.begin(), incident_net_bucket.end());
            const size_t bucket_degree = std::distance(incident_net_bucket.begin(), first_invalid_entry_it);
            const size_t tmp_incident_nets_pos = incident_nets_pos;
            incident_nets_pos += bucket_degree;
            memcpy(tmp_incident_nets.data() + tmp_incident_nets_pos,
                   incident_net_bucket.data(), sizeof(HyperedgeID) * bucket_degree);
            duplicate_incident_nets_map.clear(bucket);
          }
          
          
          // Update number of incident nets of high degree vertex
          const size_t contracted_size = incident_nets_pos - incident_nets_start;
          tmp_hypernodes[coarse_hn].setSize(contracted_size);
          
          if (deterministic) {
            // sort for determinism
            std::sort(tmp_incident_nets.begin() + incident_nets_start,
                      tmp_incident_nets.begin() + incident_nets_start + contracted_size);
          }
        }
        duplicate_incident_nets_map.free();
      }
    }

    // #################### STAGE 3 ####################
    // In the step before we aggregated hyperedges within a bucket data structure.
    // Hyperedges with the same hash/footprint are stored inside the same bucket.
    // We iterate now in parallel over each bucket and sort each bucket
    // after its hash. Parallel hyperedges are detected by comparing the pins of 
    // hyperedges with the same hash.

    // Helper function that checks if two hyperedges are parallel
    // Note, pins inside the hyperedges are sorted.
    auto check_if_hyperedges_are_parallel = [&](const HyperedgeID lhs,
                                                const HyperedgeID rhs) {
      const Hyperedge& lhs_he = tmp_hyperedges[lhs];
      const Hyperedge& rhs_he = tmp_hyperedges[rhs];
      if (isOriginalSizeUsageInParallelNetsDetectionEnabled()
          && lhs_he.originalSize() != rhs_he.originalSize()) {
        return false;
      } else if ( lhs_he.size() == rhs_he.size() ) {
        const size_t lhs_start = lhs_he.firstEntry();
        const size_t rhs_start = rhs_he.firstEntry();
        for ( size_t i = 0; i < lhs_he.size(); ++i ) {
          const size_t lhs_pos = lhs_start + i;
          const size_t rhs_pos = rhs_start + i;
          if ( tmp_incidence_array[lhs_pos] != tmp_incidence_array[rhs_pos] ) {
            return false;
          }
        }
        return true;
      } else {
        return false;
      }
    };

    for (size_t bucket = UL(0); bucket < hyperedge_hash_map.numBuckets(); ++bucket) {
      auto& hyperedge_bucket = hyperedge_hash_map.getBucket(bucket);
      std::sort(hyperedge_bucket.begin(), hyperedge_bucket.end(),
                [&](const ContractedHyperedgeInformation& lhs, const ContractedHyperedgeInformation& rhs) {
                  return std::tie(lhs.hash, lhs.size, lhs.he) < std::tie(rhs.hash, rhs.size, rhs.he);
                });

      // Parallel Hyperedge Detection
      for ( size_t i = 0; i < hyperedge_bucket.size(); ++i ) {
        ContractedHyperedgeInformation& contracted_he_lhs = hyperedge_bucket[i];
        if ( contracted_he_lhs.valid ) {
          const HyperedgeID lhs_he = contracted_he_lhs.he;
          HyperedgeWeight lhs_weight = tmp_hyperedges[lhs_he].weight();
          for ( size_t j = i + 1; j < hyperedge_bucket.size(); ++j ) {
            ContractedHyperedgeInformation& contracted_he_rhs = hyperedge_bucket[j];
            const HyperedgeID rhs_he = contracted_he_rhs.he;
            if ( contracted_he_rhs.valid &&
                 contracted_he_lhs.hash == contracted_he_rhs.hash &&
                 check_if_hyperedges_are_parallel(lhs_he, rhs_he) ) {
              // Hyperedges are parallel
              lhs_weight += tmp_hyperedges[rhs_he].weight();
              contracted_he_rhs.valid = false;
              valid_hyperedges[rhs_he] = false;
            } else if ( contracted_he_lhs.hash != contracted_he_rhs.hash  ) {
              // In case, hash of both are not equal we go to the next hyperedge
              // because we compared it with all hyperedges that had an equal hash
              break;
            }
          }
          tmp_hyperedges[lhs_he].setWeight(lhs_weight);
        }
      }
      hyperedge_hash_map.free(bucket);
    }

    // #################### STAGE 4 ####################
    // We map the global community ids to the new community ids in the coarse hypergraph.
    // This Hypergraph is then coarsened here by swapping data from temporary
    // buffers with corresponding members in this hypergraph. For the
    // incidence array, we compute a prefix sum over the hyperedge sizes in
    // the contracted hypergraph which determines the start position of the pins
    // of each net in the incidence array. Furthermore, we postprocess the incident
    // nets of each vertex by removing invalid hyperedges and remapping hyperedge ids.
    // Incident nets are also written to the incident nets array with the help of a prefix
    // sum over the node degrees.

    // map the global communities
    for ( PartitionID& community : global_communities ) {
      if ( community != kInvalidPartition ) {
        ASSERT(static_cast<size_t>(community) < communities.size());
        community = map_to_coarse_hypergraph(static_cast<HypernodeID>(community));
      }
    }

    // Setup new hypergraph _part_original_volumes, _part_sizes, _original_weighted_degrees
    for (HypernodeID id = 0; id < _num_hypernodes; ++id) {
      HypernodeID new_part_id = map_to_coarse_hypergraph(id);
      PartitionID old_part_id = _part_ids[id];
      _original_weighted_degrees[new_part_id] = _part_original_volumes[old_part_id];
    }
    // (singletone partition)
    for (HypernodeID part_id = 0; part_id < num_hypernodes; ++part_id) {
      _part_sizes[part_id] = 1;
      _part_ids[part_id] = part_id;
      _part_original_volumes[part_id] = _original_weighted_degrees[part_id];
    }
    
    // Compute number of hyperedges in coarse graph (those flagged as valid)
    std::vector<size_t> he_mapping(_num_hyperedges + 1, 0);
    for ( size_t he = 0; he < UL(_num_hyperedges); ++he ) {
      he_mapping[he + 1] = he_mapping[he] + valid_hyperedges[he];
    }
    const HyperedgeID num_hyperedges = he_mapping[_num_hyperedges];

    _hypernodes.resize(num_hypernodes);
    _weighted_degrees = std::vector<HypergraphVolume>(num_hypernodes, 0);
    _original_weighted_degrees.resize(num_hypernodes);
    _part_original_volumes.resize(num_hypernodes);
    _part_sizes.resize(num_hypernodes);
    _part_ids.resize(num_hypernodes);

    auto setup_hyperedges = [&] {
      for ( HyperedgeID id = 0; id < UL(_num_hyperedges); ++id ) {
        HyperedgeID he_mapping_value = he_mapping[id + 1] - he_mapping[id];
        if ( he_mapping_value > 0 /* hyperedge is valid */ ) {
          he_sizes[id] = tmp_hyperedges[id].size();
        } else {
          he_sizes[id] = 0;
        }
      }
      // Compute start position of each hyperedge in incidence array
      std::vector<size_t> num_pins_prefix_sum(_num_hyperedges + 1, 0);
      num_pins_prefix_sum[0] = 0;
      for ( size_t i = 0; i < UL(_num_hyperedges); ++i ) {
        num_pins_prefix_sum[i + 1] = he_sizes[i] + num_pins_prefix_sum[i];
      }
      const size_t num_pins = num_pins_prefix_sum[_num_hyperedges];
      _incidence_array.resize(num_pins);
      
      _hyperedges.resize(num_hyperedges);

      _total_volume = 0;
      // Write hyperedges from temporary buffers to incidence array
      for (HyperedgeID id = 0; id < UL(_num_hyperedges); ++id) {
        HyperedgeID he_mapping_value = he_mapping[id + 1] - he_mapping[id];
        if ( he_mapping_value > 0 /* hyperedge is valid */ ) {
          const size_t he_pos = he_mapping[id];
          const size_t incidence_array_start = num_pins_prefix_sum[id];
          Hyperedge& he = _hyperedges[he_pos];
          he = tmp_hyperedges[id];
          const size_t tmp_incidence_array_start = he.firstEntry();
          const size_t edge_size = he.size();
          std::memcpy(_incidence_array.data() + incidence_array_start,
                      tmp_incidence_array.data() + tmp_incidence_array_start,
                      sizeof(HypernodeID) * edge_size);
          he.setFirstEntry(incidence_array_start);
          _total_volume += edge_size * he.weight();
        }
      }
    };

    auto setup_hypernodes = [&] {
      // Remap hyperedge ids in temporary incident nets to hyperedge ids of the
      // coarse hypergraph and remove single-pin/parallel hyperedges.
      // (new!) At the same time compute weighted degree of each coarsened vertex
      for (HypernodeID id = ID(0); id < num_hypernodes; ++id) {
        const size_t incident_nets_start =  tmp_hypernodes[id].firstEntry();
        size_t incident_nets_end = tmp_hypernodes[id].firstInvalidEntry();
        for ( size_t pos = incident_nets_start; pos < incident_nets_end; ++pos ) {
          const HyperedgeID he = tmp_incident_nets[pos];
          HyperedgeID he_mapping_value = he_mapping[he + 1] - he_mapping[he];
          if ( he_mapping_value > 0 /* hyperedge is valid */ ) {
            tmp_incident_nets[pos] = he_mapping[he];
            _weighted_degrees[id] += tmp_hyperedges[he].weight();
          } else {
            std::swap(tmp_incident_nets[pos--], tmp_incident_nets[--incident_nets_end]);
          }
        }
        const size_t incident_nets_size = incident_nets_end - incident_nets_start;
        tmp_hypernodes[id].setSize(incident_nets_size);
        tmp_num_incident_nets[id] = incident_nets_size;
      };

      // Compute start position of the incident nets for each vertex inside
      // the coarsened incident net array
      std::vector<size_t> num_incident_nets_prefix_sum(num_hypernodes + 1, 0);
      for (size_t i = 0; i < UL(num_hypernodes); ++i) {
        num_incident_nets_prefix_sum[i + 1] = tmp_num_incident_nets[i] +
                                             num_incident_nets_prefix_sum[i];
      }
      const size_t total_degree = num_incident_nets_prefix_sum[num_hypernodes];
      _incident_nets.resize(total_degree);
      // Write incident nets from temporary buffer to incident nets array
      for (HypernodeID id = 0; id < num_hypernodes; ++id) {
        const size_t incident_nets_start = num_incident_nets_prefix_sum[id];
        Hypernode& hn = _hypernodes[id];
        hn = tmp_hypernodes[id];
        const size_t tmp_incident_nets_start = hn.firstEntry();
        std::memcpy(_incident_nets.data() + incident_nets_start,
                    tmp_incident_nets.data() + tmp_incident_nets_start,
                    sizeof(HyperedgeID) * hn.size());
        hn.setFirstEntry(incident_nets_start);
      }
    };
    
    setup_hypernodes();
    setup_hyperedges();
    _num_hypernodes = num_hypernodes;
    _num_hyperedges = num_hyperedges;
    _num_pins = _incidence_array.size();
    _total_degree = _incident_nets.size();        
    
    _k = num_hypernodes;
  }

  
  // ! Copy static hypergraph sequential
  SequentialIPHypergraph SequentialIPHypergraph::copy() const {
    SequentialIPHypergraph hypergraph;

    hypergraph._num_hypernodes = _num_hypernodes;
    hypergraph._num_removed_hypernodes = _num_removed_hypernodes;
    hypergraph._num_hyperedges = _num_hyperedges;
    hypergraph._num_removed_hyperedges = _num_removed_hyperedges;
    hypergraph._original_max_edge_size = _original_max_edge_size;
    hypergraph._num_pins = _num_pins;
    hypergraph._k = _k;
    hypergraph._total_degree = _total_degree;
    hypergraph._total_weight = _total_weight;
    hypergraph._total_volume = _total_volume;
    hypergraph._original_total_volume = _original_total_volume;
    hypergraph._disable_single_pin_nets_removal = _disable_single_pin_nets_removal;
    hypergraph._use_original_size_in_parallel_nets_detection = _use_original_size_in_parallel_nets_detection;
    hypergraph._has_original_edge_sizes = _has_original_edge_sizes;

    hypergraph._hypernodes.resize(_hypernodes.size());
    memcpy(hypergraph._hypernodes.data(), _hypernodes.data(),
           sizeof(Hypernode) * _hypernodes.size());
    hypergraph._incident_nets.resize(_incident_nets.size());
    memcpy(hypergraph._incident_nets.data(), _incident_nets.data(),
           sizeof(HyperedgeID) * _incident_nets.size());

    hypergraph._hyperedges.resize(_hyperedges.size());
    memcpy(hypergraph._hyperedges.data(), _hyperedges.data(),
           sizeof(Hyperedge) * _hyperedges.size());
    hypergraph._incidence_array.resize(_incidence_array.size());
    memcpy(hypergraph._incidence_array.data(), _incidence_array.data(),
           sizeof(HypernodeID) * _incidence_array.size());
    hypergraph._weighted_degrees.resize(_weighted_degrees.size());
    memcpy(hypergraph._weighted_degrees.data(), _weighted_degrees.data(),
           sizeof(HypergraphVolume) * _weighted_degrees.size());
    hypergraph._original_weighted_degrees.resize(_original_weighted_degrees.size());
    memcpy(hypergraph._original_weighted_degrees.data(), _original_weighted_degrees.data(),
           sizeof(HypergraphVolume) * _original_weighted_degrees.size());
    hypergraph._part_original_volumes.resize(_part_original_volumes.size());
    memcpy(hypergraph._part_original_volumes.data(), _part_original_volumes.data(),
            sizeof(HypergraphVolume) * _part_original_volumes.size());
    hypergraph._part_sizes.resize(_part_sizes.size());
    memcpy(hypergraph._part_sizes.data(), _part_sizes.data(),
            sizeof(HypernodeID) * _part_sizes.size());
    hypergraph._part_ids.resize(_part_ids.size());
    memcpy(hypergraph._part_ids.data(), _part_ids.data(),
           sizeof(PartitionID) * _part_ids.size());

    hypergraph._beta_ref = _beta_ref;
    hypergraph._gamma_ref = _gamma_ref;
    hypergraph._omega_ref = _omega_ref;

    return hypergraph;
  }

  void SequentialIPHypergraph::memoryConsumption(utils::MemoryTreeNode* parent) const {
    ASSERT(parent);
    parent->addChild("Hypernodes", sizeof(Hypernode) * _hypernodes.size());
    parent->addChild("Incident Nets", sizeof(HyperedgeID) * _incident_nets.size());
    parent->addChild("Hyperedges", sizeof(Hyperedge) * _hyperedges.size());
    parent->addChild("Incidence Array", sizeof(HypernodeID) * _incidence_array.size());
    parent->addChild("Weighted Degrees", sizeof(HypergraphVolume) * _weighted_degrees.size());
    parent->addChild("Original Weighted Degrees", sizeof(HypergraphVolume) * _original_weighted_degrees.size());
    parent->addChild("Partition Volumes", sizeof(HypergraphVolume) * _part_original_volumes.size());
    parent->addChild("Partition Sizes", sizeof(HypernodeID) * _part_sizes.size());
    parent->addChild("Partition IDs", sizeof(PartitionID) * _part_ids.size());
  }

  // ! Computes the total node weight of the hypergraph
  void SequentialIPHypergraph::computeAndSetTotalNodeWeight() {
    _total_weight = 0;
    for (HypernodeID hn = 0; hn < _num_hypernodes; ++hn) {
      if (nodeIsEnabled(hn)) {
        _total_weight += this->_hypernodes[hn].weight();
      }
    }
  }
} // namespace mt_kahypar::ds
