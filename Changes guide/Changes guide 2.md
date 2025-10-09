# Changes Guide: SoSe25

## Other changes, Ideas:
- in `computeAONParameters`, use `nodeWeightedDegree` for `ClosVol` instead of `nodeDegree` \
(as cutting edges are considered vith weights) and `totalVolume` instead of `initialTotalVertexDegree` for `vol_H` &rarr; done (for now)
- ~~in `partitioner.cpp` change `k` to the number of found communities, if `cluster && community_detection=true` (for now, I just recalculate the number of communities as done in `partitioner_output.cpp`)~~ [this way the kernel is bigger &rArr; AON IP couldn't find a good clustering and returned too many clusters...]
- in `multilevel.cpp` before uncoarsening stage check if k has changed (+correctly change it, if so) &rarr; done, not used for now
- **!!! I concider edge weights in the gain &rArr; use weighted degrees and use edge weight in _delta_cut** &rarr; ASK
- [done] enable single pin net removal (they are never cutting and volumes are tracked correctly)
- [done] set `_beta[d] = 0` if `NaN`
- [done] disable balancing for clustering in lp
- set lp to sequential (as it is bad in parallel :( [not anymore!])
- run both original sized and reduced sized versions of AON IP + a couple of random IP's [no need, as AON IP with original sizes is ~always better]
- [done] run lp-refiner in the kernel after IP [**ToDo:** Check, if  works as intended]
- [done] for edgeSizeThreshold, use `max(100, large_hyperedge_size_threshold / 10)`
- use V-cycles in multilevel (for clustering) &rarr; done

## Questions
- singleton &rarr; conductance = 1. Why?

## TODO:
- TODO: ensure, that `use_community_detection` is enabled by `aon_hypermodularity` IP [`partitioner.cpp preprocess(..)`] &rarr; done in `context.cpp sanity_check(..)`
- change `context.partition.k` in `multilevel.cpp` if `aon_hypernodularity` IP is used &rArr; done (with `fitK()` in `apply()` in `initial_partitioning_data_container.h`)
- introduce ClusteringMode to mark that `k` can be changed (?)
- run multiple AON IP in paralle (#AON = #threads) &rarr; **Problem:** parallel local search works bad: set it to sequential?
- try to snapshot and preserve original edge sizes from the start (if they are not too big: <500?)
- compare AON IP with Hypermodularity (+metrix)
- compare clustering with Hypermodularity (metric = AON &rarr; Adil has implemented the objective)

- iterate between local and global conductance in local search
- [done] merge with Adil to allow `double` gains:
    - adjust according to comments in PR #1
    - set `scaling_factor = 1` in `hypergraph_common.h`
    - [already was done when scaling was 1000] maybe than return to `delta` = sum of all local "gains" by `conductance_local` in `refineImpl()` (`label_propagation_refiner.cpp`)
- try to improve local search (gain cache? another approach?)

- **!!!** Fix a bug in IP: with `-t 8` IP fails on the asserion `partID(n) == kInvalidPartition` in `setOnlyNodeId(n, p)` &rarr; some partitions aren't empty... 

## Initial Partitioning: Hypermodularity

*(new!)*:
- `aon_hypermodularity` uses original edge sizes (if possible), no clusterpenalty
- `aon_hypermodularity_kernel` uses the edge sizes in the kernel, no clusterpenalty
- `aon_hypermodularity_bayesian` uses original edge sizes (if possible), with clusterpenalty (to lower the number of clusters) :
    + `clusterPenalty != 0` [ToDo: experiment with different (scaling) values]
    + `cluster_gain = clusterPenalty * (#clusters_before - #clusters_after)`
    - Bayesian prior regularization term from the [AON_HMLL-paper][https://arxiv.org/pdf/2101.09611]: $-#hn \cdot \log #clusters$

ToDo:
- calculate parameters ($\omega_{k0}, \omega_{k1} \beta, \gamma$) needed for the AON-Hypermodularity IP &rarr; community detection in preprocessing
- implement Hypermodularity as an IP (no parallelization in IP as multiple IP algorithms could be ran simultaneously + kernel should normally be small)
- introduce the new IP to the system:
    - register IP as an IP
    - set it as default IP for `cluster` preset (`i-runs=3`)
    - adjust the AON-parameter calculation to take place only if the IP is chosen
- try and debug

### Calculate the parameters
Reference: [commit](https://github.com/adilchhabra/mt-kahypar/commit/ab9be0777bbe77c158bf8e6f53166ea3c67ce526) My change: weighted degrees, total volume.
- `partition\partitioner.cpp`:
    - `precomputeHyperModularityParameters(&hypergraph, &context)` [my: runs only if a hypermodularity ip algo is enabled] - called by `preprocess(&hg, &context, *target_graph)` if `use_community_detection == true`
- AON-Hypermodularity block in `static_hypergraph.h, .cc`:
    - private members (copied in constructors, `contract(..)`, `copy(..)`): 
    here `_beta`, `_gamma` are the coefficients needed in the objective function
    and not the $\beta, \gamma$ from the article. \ 
    Explicitly, $\_beta[k] = \beta_k$, $\_gamma[k] = \beta_k * \gamma_k$
    ```cpp
        // AON HyperModularity Clustering Coefficients
        vec<double> _beta;                 ///< β_k
        vec<double> _gamma;                ///< β_k * γ_k
        vec<std::array<double, 2>> _omega; ///< {ω_k0, ω_k1} (ω_in, ω_out)

        
        // (15): Q(z, \Omega, d) = - \sum_k( \beta_k * (cut_k + \gamma_k * \sum_l (vol(l)^k) ) ) + J(\omega)
        //      => -\beta_k and - \beta_k * \gamma_k are the coefficient needed in the objective function
    ```
    - getters `bool hasAON()`, `double beta(d)`, `double gamma(d)`, `double omegaIn(d)`, `double omegaOut(d)`
    - procedure `void computeAONParameters(eps)`
        - uses `nodeDegree` instead of `nodeWeightedDegree` to calculate `ClusVol`
        - uses `initialTotalVertexDegree()` instead of `totalVolume()` to calculate `vol_H`
        - sets `_beta`, `_gamma`, `_omega` vectors
- mirroring interface in `dynamic_hypergraph.h`, `static_graph.h`, `dynamic_graph.h` 
  with `false` in `hasAON` and otherwise exceptions
- getters in `partitioned_hypergraph.h`, `static_hypergraph.h` and a mirroring interface in `partitioned_graph.h` and other (hyper-)graphs
    - &rArr; dummy empty double vectors in base classes to return, if AON is not supported :(

#### NaN Problem
Originally infinite values of beta and gamma were handled this way:
```cpp
_beta[d] = std::log(omega_in) - std::log(omega_out);
if (!std::isfinite(_beta[d])) {
    _beta[d] = (_beta[d] > 0.0 ? 1e3 : -1e3);
}
```
Idea was, to round infinities to +- 1000.

Problems:
1. `_beta[d]` is never negative (? &rarr; ask) &rArr; never `-Inf`
2. sometimes it's `NaN` (when one of omegas is Nan, or both are 0 or `+Inf`)

&rarr; Current solution:
0. ~~use `long double` in volume calculations to avoid `Inf` as much as possible (withput using GMP...) [Note: if `long doubles` are used for beta and gamma too, then the running time increases significantly (from <1 min to 10-20 min)]~~ [too slow, no benefit]
0. In `hypergraph_common.h`, set `AONCoefficient = double` for `_beta`, `_gamma`, `_omega` and use this alias in hypergraphs, AON IP (not yet in Objective gain computation)
1. Avoid `NaN` in omegas, by avoiding `vol_out = Inf - Inf = Nan`. If `vol_in = +Inf`, then `vol_out` should be set to `+Inf` as well \ 
    (as $(vol\_H)^d = (\sum_i vol(i))^d, vol\_in = \sum_i vol(i)^d \implies vol\_out = vol\_H^d - vol\_in$, so if `vol_in` is infinite, then `vol_out` should be at least near to `+Inf`)
2. Set `_beta[d] = 0` if it's `Nan` (when both omegas are 0, `_beta[d] = log (+Inf / +Inf) = Inf - Inf = NaN`)
```cpp
    long double vol_H_d = std::pow(vol_H, static_cast<int>(d));
    if (std::isfinite(vol_in) || std::isfinite(vol_H_d)) // 0 or Inf
      vol_out = vol_H_d - vol_in;
    else {
      ASSERT(vol_in > 0, "vol_in is not finite, but non-positive: " << vol_in);
      // 1 (Avoiding inf - inf = NaN)
      vol_out = std::numeric_limits<long double>::infinity();
    } 

    double omega_in =  static_cast<double>( (m_k[d] - cut_k[d]) / vol_in ); 
    double omega_out = static_cast<double>( cut_k[d] / vol_out );
    // (Normally) not NaN / +-Inf as m_k, cut_k < total_volume. (vol_in and vol_out != NaN, 0)

    _omega[d] = {omega_in, omega_out};
    _beta[d] = std::log(omega_in) -  std::log(omega_out);
    _gamma[d] = omega_in - omega_out;
    
    // Adjustments for Inf, Nan
    if (!std::isfinite(_beta[d])) {
      if (_beta[d] != _beta[d]) // NaN
        _beta[d] = 0; // both omegas are the same => Inf => beta = 0
      else // +Inf
        _beta[d] = (_beta[d] > 0 ? 1e3 : -1e3);
    // Idea: _beta[d] = NaN => both log omega-s are the same Inf => _beta[d] = 0
    }
```

#### Problem of huge (original) edges
If the original edge sizes are too big, the `_beta` and `_gamma` vectors are also huge. But *normally* most of their entries are 0.0 (`-1000` and `NaN` before &rarr; see *NaN Problem*).

&rarr; I track the last non-zero pair `(_beta[d], _gamma[d])` and cut all other zeroes. The getters of beta and gamma will return 0.0, if `d` is greater than the last non-zero index.

(new!) `aon_hypermodularity` IP uses the original edge sizes if available, `aon_hypermodularity_kernel` uses the edge sizes in the kernel.

#### Problem with shared between IPs local_hg 
##### Problem:
Sometimes one instance of IP gets the same local hypergraph as another IP. Because of that, partitionings could be (and sometimes are!) overwritten and even worse - merged. Happens with IPs using multiple threads (rand: in initializePertition(), AON IP: everywhere...)

To verify, I added local variable `_current_tag = -1`to `_ip_data` to store the `_tag` of the IP, which uses the local hg (analogously to the local hg itself). And made Random IP output, if `_current_tag != -1` at the start of IP or if `_current_tag != _tag` at its end.

With this modifications + `-t 10 -i-runs=100 --i-enabled-ip-algo=0 0 0 0 0 0 0 0 1 0 0 0 0` [only Rand IP], I've caught several wrong `_tag`s.

##### Reason:
1. `pool_initial_partitioner.cpp` realizes IP runs as tasks in a `tbb::task_group tg`
2. in a “global” `InitialPartitioningDataContainer ip_data` shared by all the tasks, there is a `ThreadLocalHypergraph`  (`tbb::enumerable_thread_specific<ThreadLocalHypergraph>` with a member `PartitionedHypergraph _partitioned_hypergraph`) 
3. each IP begins with getting it's thread-local hypergraph and ends by commiting the result to the global `InitialPartitioningDataContainer`
```cpp
PartitionedHypergraph& hg = _ip_data.local_partitioned_hypergraph();
...

[AON IP:] hg.initializePartition();
_ip_data.commit(/* ipName */, _rng, _tag, time);
```
4. `initializePartition()` called in parallel many things (e.g. initializes `_conductance_pq`, which uses ` vec = std::vector<T, tbb::scalable_allocator<T> >`...)
5. &rArr; an IP can be pushed from its thread into an idle state (e.g. when it waits for free threads for `parallel_for` OR for a `_pop_lock` in `commit(..)`, if `--deterministic=true` was set...)
6. &rArr; later IP potentially resumes on some other thread (or a new IP starts on a thread of an idle IP)
7. &rArr; two IPs have references to the same `_local_hg`
8. &rArr; their partitionings can be merged / rewritten :(

##### Current solution
- In `pool_initial_partitioner.cpp` instead of `tbb::task_group tg`, use `tbb::parallel_for` to run IPs in parallel (as `tbb::task_group` can switch tasks between threads, but `tbb::parallel_invoke` doesn't) [Potentially overkill, none the less it works]
- Use a `tbb::task_arena` limited to 1 thread, to run IP on a single thread:

```cpp
 if ( _ip_data.should_initial_partitioner_run(_ipName) ) {
    /* Concurrency Problem:
     * - IP gets a local to its thread hypergraph
     * - if an IP used multiple threads, and there are no more threads available,
     *   tbb can make this IP wait and start another IP in the same thread
     *   (the same happens with blocked by a mutex IPs)
     * => two IPs can share the same local hypergraph
     * => they can override / merge each other's partitioning
     * 
     * Current solution: limit IP to one thread
    */
    tbb::task_arena limited_arena(1);
    limited_arena.execute([&]{
        [AON IP]
    });
 }
``` 
##### Another tried solution [Reverted]
1. I've defined `SequentialIPHypergraph`:
- sceletone: `StaticHypergraph` (`SequentialIPHypergraph` is a `friend` of `StaticHypergraph`)
- differences:
    - no multithreading, `vec`, parallel initializations of `Array`s
    - `_part_ids` instead of `_community_ids`
    - `_part_sizes`, `_part_original_volumes`
    - initializes from `StaticHypergraph` with a singleton partitioning
    - in `contract(global_communities)`:
        - contracts all clusters in the hypergraph itself &rarr; no returned hypergraph
        - translates values in `global_communities` to the new IDs

2. For that, an additional class was needed: `SequentialBuchetMap` (analog to `ConcurrentBuchetMap`).

3. I've rewritten AON IP to use only `SequentialIPHypergraph`. As it doesn't have `_pin_count_in_part` ($O(k \cdot #he) = O(#hn \cdot #he)$...), I've used the Hypermodularity.jl approach with `notsame`.
4. Also, to remove concurrenct form the IP, I've changed `conductance_pq.h` to use `priority_queue_sec.h` - an analogon of `priority_queue.h`, but with `std::vector` instead of `vec`.

Sadly, this only reduced the frequency of wrong tags, but didb't solve the problem &rArr; **these changes are reverted**

**Saved files:**
- in `mt-kahypar/partition/initial_partitioning`
    - `aon_ip_temp.h, .cc`
- in `mt-kahypar/datastructures`
    - `priority_queue_seq.h`
    - `sequential_ip_hypergraph.h, .cpp`
    - `sequential_bucket_map.h`
    - `conductance_pq_temp.h`
**Reverted:** 
- `CmakeLists.txt` (added the new `.cpp` files)
- `StaticHypergraph`:
```cpp
class SequentialIPHypergraph;
class StaticHypergraph { ... friend SequentialIPHypergraph; ... }
```

### Implement AON-Hypermodularity IP
Original Algorithm: [Generative hypergraph clustering: from blockmodels to modularity](https://arxiv.org/pdf/2101.09611)

#### Code Structure 
("Template": `mt-kahypar/partition/initial_partitioning/random_initial_partitioner.h, .cpp`)

`aon_hypermodularity_initial_partitioner.h, .cpp`:
- Interface for customized calling:
```cpp 
  void partitionImpl(const HypernodeID edgeSizeThreshold, 
                     const long long maxNumIter, 
                     const double eps, 
                     const double clusterPenalty,
                     const bool randomize,
                     bool useOriginalEdgeSizes,
                     const InitialPartitioningAlgorithm ip_name);
  
  void partitionImpl() final {
    partitionImpl(
      1e3 /* edgeSizeThreshold */,
      1e2 /* maxNumIter */,
      1e-8 /* eps */,
      0.0 /* clusterPenalty */,
      true /* randomize */,
      true /* useOriginalEdgeSizes */,
      InitialPartitioningAlgorithm::aon_hypermodularity /* ip_name */
    );
  }
```
- Main code:
```cpp
      if (! useOriginalEdgeSizes || ! H.hasOriginalEdgeSizes()) {
        useOriginalEdgeSizes = false;
        if ( ! H.hasOriginalEdgeSizes() )
          LOG << "No snapshot of original edge size found";
        // save current edge sizes, weighted degrees and total volume
        H.snapshotOriginalEdgeSizes();
        H.snapshotOriginalWeightedDegreesAndTotalVolume();
        H.useOriginalSizeInParallelNetsDetection(true); // otherwise gain is incorrect
      }
      if (_context.partition.verbose_output) {
        if (useOriginalEdgeSizes)
          LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: Using original edge sizes: # he " << H_new.initialNumEdges() << "# hn " << H_new.initialNumNodes();
        else
          LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: Using current edge sizes. # he " << H_new.initialNumEdges() << "# hn " << H_new.initialNumNodes();
      }
      H.enableSinglePinNetsRemoval(); // single pin nets are never cutting

      //          1. Singleton initial partitioning
      //                         <...>
      //          2. AllOrNothingHMLL: Louvain Cycle   
      double total_gain = 0.0;
      long long counter = 1;
      while (z_changed) {
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

        /** ------------------ Louvain Step: ------------------
         *  - Nodes are moved to neighbouring partitions as
         *    long as it improves the modularity gain;
         *  - map_z is updated accordingly.
         */
        AONCoefficient new_gain = louvainStep(H_new, H_new_partitioned, map_z, clusterSizes,
                                beta, gamma, edgeSizeThreshold, maxNumIter, eps, clusterPenalty, randomize);
        if (_context.partition.verbose_output)
          LOG << "AON IP [" << V(_ipName) << V(_tag) << "]: step " << counter << " gain = " << new_gain;
        total_gain += new_gain;
        z_changed = (new_gain > eps);

        /** --------------------- Expand: ---------------------
         *  - If H_new_partitioned is still in a singleton 
         *    partition, false is returned;
         *  - otherwise, the community structure on H_new is
         *    updated according to the new partition on 
         *    H_new_partitioned;
         *  - z (communities for H) is updated via map_z;
         *  - true is returned.
         */
        
        expand(H, H_new, H_new_partitioned, map_z, z);
      } while (z_changed && counter++ < 10 && H_new_partitioned.k() > 2); // max. 10 iterations

      if (_context.partition.verbose_output)
        LOG << "AON IP [" << V(_ipName) << V(_tag) << "] finished Louvain with total gain " << total_gain 
            << " and " << H_new_partitioned.k() << " clusters.";

      //             3. Finalize Partitioning
      //                        <...>

```

`aon_hypermodularity_kernel_initial_partitioner.h, .cpp`:
- two attributes:
    - `AONHypermodularityInitialPartitioner _aon_ip`
    - `const Context& _context` (to get `k_context.partition.large_hyperedge_size_threshold` for the `edgeSizeThreshold` parameter)
- calls `aon_hypermodularity` with `useOriginalEdgeSizes=false`


`aon_hypermodularity_bayesian_initial_partitioner.h, .cpp`:
- two attributes:
    - `AONHypermodularityInitialPartitioner _aon_ip`
    - `const Context& _context` - to get:
        + `_context.partition.large_hyperedge_size_threshold` for the `edgeSizeThreshold` parameter
        + `_context.partition.initial_num_nodes` for the `clusterPenalty` parameter
- calls `aon_hypermodularity` with `clusterPenalty = initial_num_nodes / 100`

#### Needed Additional Functionality (Static Hypergraph, Partitioned Hypergraph, Context)

1. [done] maintain original edge size in Static Hypergraph + **mirroring public interface in partitioned** / static / dynamic (hyper)graphs, **in both `copy(..)` in `static_hypergraph.cpp` copy `_original_max_edge_size`and in the constructor**:
```cpp
   class Hyperedge { // copy constructor in contract() should preserve all information -> no changes needed
    /// ...
    private:
        size_t original_size; // Number of pins at the moment of snapshot
        size_t originalSize() const;
        void setOriginalSize(size_t size);
    };

    class StaticHypergraph {
        // ...
        public:
            bool hasOriginalEdgeSizes() const;
            void snapshotOriginalEdgeSizes(); // saves sizes, manually computes _original_max_edge_size (as _max_edge_size could not be initialized yet)
            HypernodeID originalEdgeSize(HyperedgeID e) const;
            HypernodeID originalMaxEdgeSize(HyperedgeID e) const; 
        private:
            HypernodeID _original_max_edge_size; // set at the moment of snapshot
    };
```
2. coarsest_underlying_hg should forget initial weighted degrees etc (as otherwise the original edge sizes are too big => delta to slow) /
    ~~-> use factory to copy the hypergraph (check before doing!) [doesn't work]~~
    + add forgetting functions to `static_hypergraph.h`, `partitioned_hypergraph.h` and mirroring interface to others:
    ```cpp
    class StaticHypergraph, DynamicHypergraph {
    public:
        // together for their consistency
        void snapshotOriginalWeightedDegreesAndTotalVolume();
    private:
        void snapshotOriginalWeightedDegrees();
        void snapshotOriginalTotalVolume();
    };

    class PartitionedHypergraph, PartitionedGraph {
    public:
        // together for their consistency
        // recomputes conductance_pq if needed and uses original stats
        void snapshotOriginalWeightedDegreesAndVolumes();
    private:
        void snapshotOriginalWeightedDegrees();
        void snapshotOriginalTotalVolume();
        void snapshotOriginalPartVolumes();
    };
    ```
3. `enableSinglePinNetsRemoval()` in `StaticHypergraph` (+ mirroring in dynamic + graphs) - to remove all single-pin nets in the coarsest hypergraph \ 
    ~~**Rationale**: ~~\
    ~~The edge size of the original hypergraph is potentially too big to compute Hypermodularity gains~~ \
    &rArr; ~~I run AON-Hypermodularity on the actual stats of the coarsest hypergraph~~ \
    &rArr; ~~single-pin nets are never cut~~ \
    &rArr; ~~are interesting only for value~~ \
    &rarr; ~~can be removed after saving original volumes and weighted degrees~~ \
    [actually irrelevant, as parallel single pin nets are removed]
4. `useOriginalSizeInParallelNetsDetection(bool yes)` in `StaticHypergraph` (+ mirroring in dynamic + graphs) - to stop removal of parallel nets of different original sizes (otherwise the gain is incorrect) (**in both `copy(..)` in `static_hypergraph.cpp` copy `_use_original_size_in_parallel_nets_detection`and in the constructor**):
    ```cpp
    void useOriginalSizeInParallelNetsDetection(bool yes = true) {
        _use_original_size_in_parallel_nets_detection = yes;
    }
    bool isOriginalSizeUsageInParallelNetsDetectionEnabled() const {
        return _use_original_size_in_parallel_nets_detection;
    }

    ```
5. ~~Saving k to update (decrease) later, as AON IP decreases k (sometimes significantly):~~ [doesen't work as Pool IP discards the hg info and saves only partition &rArr; saved k is lost] \
    Fitting k to the minimal possible number without changing the partIDs:
    - Rationale: before IP k is set to `num_hypernodes`, IP can reduce it greatly (And AON IP uses the smallest possible PartID's!!!)
    - Implementation: find the maximal used part ID and set k to it
    - if `cluster` mode is used , called in `apply()` before `_partitioned_hg.initializePartition()` for the best IP result (`initial_partitioning_data_container.h`)

    - `partitioned_hypergraph.h, partitioned_graph.h`:
        + \+ `fitK()`
```cpp

  // ! Fits k before calling initializePartition()
  void fitK() {
    // accumulate in parallel the maximal used part ID
    PartitionID maxUsedPartID = tbb::parallel_reduce(/* max */);
    PartitionID actualK = 1 + maxUsedPartID;
    ASSERT(actualK <= _k);
    if (actualK < 2) {
      actualK = 2;
      LOG << RED << "PartitionedHypergraph::fitK() - Warning: only one cluster found: "
             "actualK = " << actualK << ", setting it to 2" << END;
    }
    if (_k != actualK) {
      setK(actualK);
      LOG << "PartitionedHypergraph::fitK() - Fitted k to " << actualK;
    }
  }

```
6. [done] Save initial number of nodes in `context.partition.initial_num_nodes` (to use it as `clusterPenalty` in `aon_hypermodularity_bayesian` IP) - analogously to `context.partition.initial_k`:
    - `context.h`: 
        + `HypernodeID initial_num_nodes = -1;` to `struct PartitionContext`
    - `partitioner.cpp`:
        + in `setupContext(...)` set `context.partition.initial_num_nodes = hypergraph.initialNumNodes();`
    - `sql_plottools_serializer.cpp`:
        + in `void SQLPlottoolsSerializer::writePartitionContext(...)` print out `initial_num_nodes`


#### Implementation Details
1. Louvain `Collapse(..)` and `Expand(..)`:
    ![Algorithm 3](<Algorithm 3: AllOrNothingHMLL.png>)
    - **Main idea**: create a partitioned hg from singletons = collapse
    - in `collapse(..)` save mapping to `map_z`: \
    	`map_z[communityID(collapsed_hn)] = collapsed_hn`, initialize `clusterSizes` with ones [for cluster penalty]
    - in `expand(..)` adjust `z` (*= the best found partitioning of the given coarsest hypergraph*), contract static hg [here, as in `collapse` it's unnecessary at the first time and potentially takes too long]: \
	    ```
        for hn in H:
            z_new[hn] = map_z[z[hn]]
        z = z_new
        ```
2. Louvain step:
   - before each move update `map_z`:
	`map_z[community_id[hn]] = <new_partition_id>` \
	~~`map_z[community_id[hn]] = map_z[community_id[hn_of_new_label]]`~~
	~~(`hn_of_new_label` should be equal to the new `CommunityID A`)~~
    ![Algorithm 4](<Algorithm 4: AONLouvainStep.png>)
    \+ `eps = 1e-8`, `maxNumIter = 1000`, `randomize = true` - default parameters:
    - `if best_gain > eps` the move is made. Otherwise never stops: 
        >    ... \
        >    Louvain: round 618 \
        >    Louvain: node 0 \
        >    Louvain: node 1000 \
        >    Louvain: node 2000 \
        >    Louvain: node 2874 -> 2874 (gain: **2.47727e-46**) \
        >    Louvain: node 3000 \
        >    Louvain: node 3853 -> 4414 (gain: **2.47727e-46**) \
        >    ...
    - `if randomize`: the nodes are contracted in a random order:
        ```cpp
        double total_gain = 0.0;
        if (randomize) {
          //              shuffle nodes
          while (improving && (iter++ < maxNumIter)) {
            improving = false;
            double gain = 0.0;
            for (const HypernodeID &i : nodes) {
              gain += louvainStepForANode(i, ...);
            }
            total_gain += gain;
            improving = (gain > eps);
          } 
        } else {
          // the same, but without shuffling
        }
        ```

3. Q_AON Gain:
    ![Algorithm 5](<Algorithm 5: QAON gain.png>)
    Here:
    - $d$ - the original (by me weighted) degree in $H$
    - $vol$ - the original (by me weighted) volume in $H$
    - $\beta_k$ is stored in `_beta`
    - $\beta_k \cdot \gamma_k$ is stored in `_gamma`
    - $k^{\_}$ is the maximal edge size in `H`
    - $s^{\_}$ is the original edge size in $H$ *(removal of parallel edges is not a problem)*
    
    **!!! I concider edge weights in the gain &rArr; use weighted degrees and use edge weight in _delta_cut** &rarr; ASK

    [For now] I introduced `const HypernodeID edgeSizeThreshold = std::max<HypernodeID>(_context.partition.large_hyperedge_size_threshold / 10, 100) /* edgeSizeThreshold */,` parameter to avoid too big loops in gain computation and (more importantly) calculating with infinities (as inf - Inf = NaN &rarr; the gain is `NaN` &rArr; move isn't done)

    Cluster penalty is added to the `QAONGain`. `clusterPenalty > 0` intensifies the reduction of the number of clusters.

1000. [old] ~~The underlying hypergraph `H` can have too many single-pin hyperedges &rArr; I contract it's singleton communities and after that disabled single-pin nets removal~~ [not a problem, as single pin nets are removed now]

### Introduce of the new IPs to the framework

1. Adjust the parameters of `cluster` preset:
    - set `k` to 32 instead of 2 (Adil has done so -> TODO: ask):
        - `mt-kahypar\partition\context.cpp` in `Context::setupContractionLimit(total_hypergraph_weight)`: \
        for `cluster` preset, set `coarsening.contraction_limit` to `coarsening.contraction_limit_multiplier * 32` instead of `coarsening.contraction_limit_multiplier * partition.k`.
        - `mt-kahypar\partition\partitioner.cpp` in `setupContext(& hypergraph, & context, *target_graph)`: \ 
        set `k=32` for `cluster` preset
    - [old changes guide] analog. to `singleton` introduce `aon_hypermodularity`, `aon_hypermodularity_kernel`, `aon_hypermodularity_bayesian` IPs:
        - in `config/`:
            - `cluster_preset.ini`: 
                - set `i-enabled-ip-algos=1` only for `aon_hypermodularity`, `aon_hypermodularity_kernel`, `aon_hypermodularity_bayesian`
            - `large_k_preset.ini`: 
                \+ `i-enabled-ip-algos=0` for `aon_hypermodularity`, `aon_hypermodularity_kernel`, `aon_hypermodularity_bayesian`
            - (all other `.ini` use `initial_partitioning: i-mode=rb` [recursive bipartitioning] &rArr; no changes)
        - `mt-kahypar/`:
            - in `partition/`:
                - `context_enum_classes.h`:
                    - in `InitialPartitioningAlgorithm`:
                    ```cpp
                    enum class InitialPartitioningAlgorithm : uint8_t {
                    ...
                    aon_hypermodularity = 10,
                    aon_hypermodularity_kernel = 11,
                    aon_hypermodularity_bayesian = 12,
                    UNDEFINED = 13
                    };
                    ```
                - `context_enum_classes.cpp`:
                    - in `operator<<(or, algo)` and `initialPartitioningAlgorithmFromString(algo)` add transmations string <-> `InitialPartitioningAlgorithm` for `aon_hypermodularity`, `aon_hypermodularity_kernel`
                    - in `isHypermodularityIP(algo)` return `true` for `aon_hypermodularity`, `aon_hypermodularity_kernel`, `aon_hypermodularity_bayesian`
                - in `initial_partitioning/`: 
                    - \+ `aon_hypermodularity_initial_partitioner.h, .cpp`
                    - \+ `aon_hypermodularity_kernel_initial_partitioner.h, .cpp` 
                    - \+ `aon_hypermodularity_bayesian_initial_partitioner.h, .cpp`
                    - **!!!** `CMakeLists.txt`: \+ `aon_hypermodularity_initial_partitioner.cpp`, `aon_hypermodularity_kernel_initial_partitioner.cpp`, `aon_hypermodularity_bayesian_initial_partitioner.cpp`
                - in `registries/`:
                    - `register_initial_partitioning_algorithms.cpp`
                    - \+ `#include "../initial_partitioning/aon_hypermodularity_initial_partitioner.h"`, `#include "../initial_partitioning/aon_hypermodularity_kernel_initial_partitioner.h"`, `#include "../initial_partitioning/aon_hypermodularity_bayesian_initial_partitioner.h"`
                    - \+ define `AONHypermodularityPartitionerDispatcher`, `AONHypermodularityKernelPartitionerDispatcher`, `AONHypermodularityBayesianPartitionerDispatcher`
                    - in `register_initial_partitioning_algorithms()`: \+ register `AONHypermodularityPartitionerDispatcher`, `AONHypermodularityKernelPartitionerDispatcher`, `AONHypermodularityBayesianPartitionerDispatcher`
            - in `io/`:
                - `command_line_options.cpp`: 
                    - by `"i-enabled-ip-algos"` example add `aon_hypermodularity`, `aon_hypermodularity_kernel`, `aon_hypermodularity_bayesian` IP (and change the number of IP-algos at the end of the example)
                - `presets.cpp`:
                    - in `load_large_k_preset()`: by`// main -> initial_partitioning` add entry for `aon_hypermodularity` (`"0"`), `aon_hypermodularity_kernel` (`"0"`), `aon_hypermodularity_bayesian` (`"0"`)
                    - in `load_clustering_preset()`: 
                        ```cpp
                         "0" // singleton" IP
                         "1" // aon_hypermodularity
                         "1" // aon_hypermodularity_kernel (always worse?)
                         "1" // aon_hypermodularity_bayesian
                        ```
    - `config/cluster_preset.ini` and `mt-kahypar/io/presets.cpp`: in `# main -> initial_partitioning` set `i-runs=2` istead of 10 for `cluster` (as AON-hypermodularity is ~~deterministic~~ randomized and runs `i-runs * t` times)

3. `sanityCheck(*target_graph)` in `context.cpp`:
    - ~~adjust~~ conductance checks to allow hypermodularity IPs (via `isHypermodularityIP(ip)`)
    - ensure (via `usesHypermodularityIP(ip)`), that `use_community_detection` is enabled if a hypermodularity IP is used

4. ~~[Idea] Change `context.partition.k` if it changed after IP (due to `aon_hypernodularity`)  &rarr; not done, as `new_k` shouldn't be greater than the number of nodes~~
4. change `context.partition.k` in `multilevel.cpp` if it `aon_hypernodularity` IP is used \ 
[analog. to `cluster` + `singleton` &rArr; `new_k = #_nodes`] 

#### Define helper functions to detect hypermodularity IP

`mt-kahypar/partition/context_enum_classes.h, cpp`:
- \+ `bool isHypermodularityIP(InitialPartitioningAlgorithm algo)`

`mt-kahypar/partition/context.h`:
- \+ `bool usesHypermodularityIP() const`

### Run Hypermodularity IPs in parallel
Hypermodularity IP heavily depends on the visiting order of the nodes.
&rArr; we run multiple hypermodularity IPs in parallel.

**For now for hypermodularity IPs, I multiply the number of runs by the number of threads.**
`--i-runs=2`

`mt-kahypar/partition/initial_partitioning/pool_initial_partitioner.cpp`:
```cpp
template<typename TypeTraits>
void Pool<TypeTraits>::bipartition(...) {
  //                        <..>
  // Push the runs of the different initial partitioning algorithms into a task list
  for ( uint8_t i = 0; i < static_cast<uint8_t>(InitialPartitioningAlgorithm::UNDEFINED); ++i ) {
    if ( context.initial_partitioning.enabled_ip_algos[i] ) {
      auto algorithm = static_cast<InitialPartitioningAlgorithm>(i);
      size_t runs = context.initial_partitioning.runs;
      if (isHypermodularityIP(algorithm)) {
        runs *= context.shared_memory.original_num_threads;
        LOG << "Running " << static_cast<InitialPartitioningAlgorithm>(algorithm) 
            << " IP " << runs << " time(s) as a Hypermodularity IP";
      }
      for ( size_t j = 0; j < runs; ++j ) {
        _ip_task_lists.emplace_back(algorithm, rng(), tag++);
      }
    }
  }
  //                        <...>
}
```

### Problems

## AON Hypermodularity Objective (for comparimg)
Reference: Adil's commit [ab9be07](https://github.com/adilchhabra/mt-kahypar/commit/ab9be0777bbe77c158bf8e6f53166ea3c67ce526#diff-3a1538bb7d62a2dd66446fc81ad9973acdfafbc028c374139a2e1977a98b89b3)


### Registration
- `include/`:
    - `mtkahypar.h` [instead of `libmtkahypar.h`]:
        + \+ `mt_kahypar_aon_hypermodularity(phg)`
    - `mtkahypartypes.h` [instead of `libmtkahypartypes.h`]:
        + \+ `mt_kahypar_objective_t AON_HYPERMODULARITY`
    - `lib_generic_impls.h` [my]:
        + \+ `aon_hypermodularity(phg)`
- `lib/mtkahypar.cpp` [instead of `lib/libmtkahypar.cpp`]:
    + \+ in `mt_kahypar_set_context_parameter(..)`, `mt_kahypar_set_partitioning_parameters(..)`, `mt_kahypar_get_objective(..)`,  add a case for `aon_hypermodularity`
    + \+ `mt_kahypar_aon_hypermodularity(phg)` - returns value
-  `mt-kahypar/`:
    - `io/`:
        - `sql_plottools_serializer.cpp`:
            - print `"aon_hypermodularity="` and `metric::quality`
        - `partitioning_output.cpp`:
            - print in `printObjectives(..)`
    - `partition/`:
        - `context.cpp`:
            - in `setupGainPolicy()` add case for `aon_hypermodularity` objective
        - `context_enum_classes`:
            - `.cpp`: add (2) cases for `aon_hypermodularity` objective 
            - `.h`: add `aon_hypermodularity` to 2 enums 
        - `metrics.cpp`:
            + \+ `struct ObjectiveFunction` for `aon_hypermodularity`
            - \+ `Gain aonVolumeTerm(phg)` - a help function (**ToDo: change to Gain = double**)
            - \+ add `aon_hypermodularity` case to `Gain compute_objective_parallel(phg)`, `Gain compute_objective_sequentially(phg)`, `quality(..)`, `contribution(..)`
        - `refinement/`: 
            - `gains`:
                - \+ `aon_hypermodularity/`:
                    + \+ `aon_hypermodularity_attributed_gains.h`
                    + \+ `aon_hypermodularity_gain_computation.h`
                - `bipartitioning_policy.h`:
                    - `useCutNetSplitting = true`, `nonCutEdgeMultiplier = 1` for `aon_hypermodularity`
                - `gain_cache_ptr.h`:
                    - add `aon_hypermodularity` case in `applyWithConcreteGainCache(..)`, `applyWithConcreteGainCacheForHG(..)`, `constructGainCache(..)` (`Km1GainCache` **!!**)
                - `gain_definitions.h`:
                    - include `aon_hypermodularity` (2)
                    case in `Gain::compute(..)`
                    - \+ `AONHypermodularityGainTypes`:
                        ```cpp
                        struct AONHypermodularityGainTypes : public kahypar::meta::PolicyBase {
                            using GainComputation = AONHypermodularityGainComputation;
                            using AttributedGains = AONHypermodularityAttributedGains;
                            using GainCache = Km1GainCache;
                            using DeltaGainCache = DeltaKm1GainCache;
                            using Rollback = Km1Rollback;
                            using FlowNetworkConstruction = Km1FlowNetworkConstruction;
                        };
                        ```
                    - add `aon_hypermodularity` to `GainTypes`, `_LIST_HYPERGRAPH_COMBINATIONS`, `_INSTANTIATE_CLASS_MACRO_FOR_HYPERGRAPH_COMBINATIONS`, `SWITCH_HYPERGRAPH_GAIN_TYPES`
        - `registries/register_policy.cpp`:
            - register `aon_hypermodularity` under `Gain Type Policies`

- [my] `python/module.cpp`: 
    - \+ add new objective in `### Enum Types ###`
    - \+ register `aon_hypermodularity` with calculating function

### Implementation


- Save betas and gammas:
    -  `mt-kahypar/datastructures/hypergraph_common.h`: add `const *beta_vec, *gamma_vec` and `max_edge_size` to `SynchronizedEdgeUpdate`:
    ```cpp
    struct SynchronizedEdgeUpdate {
      //                  <...>
      // (new) For AON-Hypermodularity:
      //! Pointers / views to the *global* βₖ and γₖ arrays of the **current
      //! hierarchy level**.  They never change during the whole node move, so a
      //! (const) pointer is enough and avoids copying the vectors for every edge.
      const vec<double>* beta_vec        = nullptr;    // same object for all edges
      const vec<double>* gamma_vec       = nullptr;
      HyperedgeID original_edge_size;
      HyperedgeID hn_degree;
      HypergraphVolume original_volume_from_after;
      HypergraphVolume original_volume_to_after;
      HypergraphVolume original_weighted_degree;
      HyperedgeID max_edge_size = kInvalidHyperedge;
    };
    
    ```
    - set in `partitioned_hypergraph.h`, `changeNodePart(..)`

- snapshot original edge sizes for AON Hypermodularity: in `partitioner.cpp`, `precomputeHyperModularityParameters()`
- disable error of missing `-k / --blocks`, `-e / --epsilon` parameter, if `clustering` is used 
    - `mt-kahypar/io/command_line_options.cpp`: `processCommandLineInput(...)`
    - `lib/mtkahypar.cpp`: `mt_kahypar_set_partitioning_parameters(...)`:
    ```cpp
      c.partition.k = num_blocks > 1 ? num_blocks : 32;
      c.partition.epsilon = epsilon >= 0.0 ? epsilon : std::numeric_limits<double>::max();
    ```
- help functions in `gain_computation_base.h`:
    + \+ `aonVolume(sync_update)` calling `AttributedGains::volumeDelta()`
    + \+ `setLocalDelta(moveGain)` (**ToDo: check, where / if used**)
    ```cpp
    void setLocalDelta(Gain moveGain) {
      _deltas.local() += moveGain;
    }
    ```
- to print the resulting double AON Hypermodularity: 
    - `long double compute_double_aon_hypermodularity(phg)` in `metrics.h, cpp` [uses `long double` to caculate volumes to avoid `Inf` as much as possible]
    - \+ `define AON_HYPERMODULARITY_DOUBLE` in `metrics.cpp`
    - print in `printObjectives` in `partitioning_output.cpp`


## Configuration of cluster preset
Use lp-refiner on the kernel after IP to improve conductance of the found clusters. *[Note: potential problem: it's almost always beneficial for `conductance_global` to merge clusters &rArr; lp refiner can merge too many clusters together]*

Bug fix: `setNodePart(..)` was used in `performRefinementOnPartition(..)` &rArr; I call `needsConductancePriorityQueue()` afterwards.


Set in ``config/cluster_preset.ini`` (and analog. in `mt-kahypar/io/presets.cpp`):
```ini
# main -> initial_partitioning
i-runs=2
i-enabled-ip-algos=1 # aon_hypermodularity
i-perform-refinement-on-best-partitions=true
i-fm-refinement-rounds=0
i-lp-maximum-iterations=0
i-lp-initial-block-size=0
...
# main -> initial_partitioning -> refinement -> label_propagation
i-r-lp-type=label_propagation
i-r-lp-maximum-iterations=5
...
```
**ToDo:** Check, if instead of this, `i-lp-maximum-iterations=0` with no `ip -> refinement -> lp` works as I intended

## Label Propagation

### Problems:
- more than 3 threads &rArr; bad results :(
    - **ToDo (?)**: set a sequential flag for LP to be able to run multiple IP in parallel
- [solved for now] lp refiner summs up all Attributed gains of the moves in a `delta`. With `HyperedgeWeight` gains I get an overflow &rArr; `refineImpl(..)` returns `detla > 0  = false` and the refinement round finishes:
    - ~~Current solution: recalculate `delta` with the actual conductance gain (To be changed (?) after introducting `double` gains)~~ [was already changed after scaling * 1000]
- [solved for now] didn't move the last node from it's cluster what resulted in a bad conductance:
    - Current solution: practically disable `isBalanced(..)` (`metrics.cpp`) check called in `LabelPropagationRefiner<GraphAndGainTypes>::labelPropagationRound(..)` (it wasn't balanced &rArr; called rebalancer &rArr; sometimes / mostly (?) finished LP) 
    ```cpp
    template<typename PartitionedHypergraph>
    bool isBalanced(const PartitionedHypergraph& phg, const Context& context) {
      size_t num_empty_parts = 0;
      bool acceptable_part_weights = true;
      for (PartitionID i = 0; i < context.partition.k; ++i) {
        if (phg.partWeight(i) > context.partition.max_part_weights[i]) {
          acceptable_part_weights = false;
          break;
        }
        if (phg.partWeight(i) == 0)
          num_empty_parts++;
      }
      return (acceptable_part_weights && context.partition.preset_type == PresetType::large_k) ||
             (acceptable_part_weights && num_empty_parts <= phg.numRemovedHypernodes()) ||
             (context.partition.preset_type == PresetType::cluster 
              && (phg.k() - num_empty_parts > 1));
    ```
- no gain cache for conductance &rarr; we use `CutGainCache`, which allocates an array of the size `original_num_nodes * k` &rArr; `segfault` on big hypergraphs (`circuit5M`):
    - **ToDo**: implement a (dummy?) gain cache [done]
- if only one cluster is left, modes aren't moved out as they are not border nodes \
 &rarr; changed `LabelPropagationRefiner<GraphAndGainTypes>::moveVertex(..)` to allow moving non-border nodes if `clustering` os used

## Side trip: GMP
We planned to use GMP to avoid `Inf` and `NaN` problems in AON Hypermodularity. But for now, we decided to use `double` and avoid too long running times (as even `long double` increased running time drastically and `GMP` seems to be even slower).

### Build: CMake
As `GMP` is not installed on the system, I needed to build it from source. \ 
*[The system didn't have many tools (at least updated) &rArr; the build process is really messy...]*

```cmake
# dependencies
...
option(KAHYPAR_DOWNLOAD_GMP "Download GMP automatically." OFF)

...

#################################################################
## Setup of dependencies                                       ##
#################################################################
...
set(KAHYPAR_GMP_VERSION          v6.2.1)

...
### [Old code - here for reference only]
if(KAHYPAR_DOWNLOAD_TBB)
  # Download TBB library
  ...
endif()
target_link_libraries(MtKaHyPar-Include INTERFACE TBB::tbb TBB::tbbmalloc)
### [End of old code - here for placing reference only]

# Download GMP Library
if(KAHYPAR_DOWNLOAD_GMP)
  # Set the GMP install directory
  set(GMP_INSTALL_DIR "${CMAKE_BINARY_DIR}/gmp")

  # Download GMP source code
  message(STATUS "Downloading GMP Source: version 6.2.1")
  FetchContent_Declare(
      GMP EXCLUDE_FROM_ALL SYSTEM
      URL https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
      DOWNLOAD_COMMAND wget --no-check-certificate https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz -O ${CMAKE_BINARY_DIR}/_deps/gmp-6.2.1.tar.lz
  )

  # Handle the unpacking, configuring, and building manually
  FetchContent_GetProperties(GMP)
  if(NOT GMP_POPULATED)
      FetchContent_Populate(GMP)

      set(GMP_SOURCE_DIR "${CMAKE_BINARY_DIR}/_deps")
      message(STATUS "Extracting GMP Source")
      # Unpack the tarball manually
      execute_process(
          COMMAND ${CMAKE_COMMAND} -E tar xvf ${GMP_SOURCE_DIR}/gmp-6.2.1.tar.lz
          WORKING_DIRECTORY ${GMP_SOURCE_DIR}
          OUTPUT_QUIET
      )

      message(STATUS "Configuring GMP")
      # Run configure for GMP
      execute_process(
          COMMAND ./configure --enable-cxx --prefix=${GMP_INSTALL_DIR}
          WORKING_DIRECTORY ${GMP_SOURCE_DIR}/gmp-6.2.1
          OUTPUT_QUIET
      )

      message(STATUS "Building GMP")
      # Build GMP
      execute_process(
          COMMAND make -j8
          WORKING_DIRECTORY ${GMP_SOURCE_DIR}/gmp-6.2.1
          OUTPUT_QUIET
      )

      message(STATUS "Installing GMP")
      # Install GMP
      execute_process(
          COMMAND make install
          WORKING_DIRECTORY ${GMP_SOURCE_DIR}/gmp-6.2.1
          OUTPUT_QUIET
      )
      message(STATUS "GMP Installed Successfully at ${GMP_SOURCE_DIR}/gmp-6.2.1")
  endif()

  # Set the include and library directories for GMP
  set(GMP_INCLUDE_DIR "${GMP_INSTALL_DIR}/include")
  set(GMP_LIBRARIES "${GMP_INSTALL_DIR}/lib/libgmp.a")

  message(STATUS "GMP Include: ${GMP_INCLUDE_DIR}, GMP Library: ${GMP_LIBRARIES}")

else()
#  # Find system GMP library - !!! didn't seem to work properly !!!
#  find_package(GMP REQUIRED)
#  if(NOT GMP_FOUND)
#    message(FATAL_ERROR "
#      GMP library not found or current GMP version is too old. Install GMP on your system
#      or add -DKAHYPAR_DOWNLOAD_GMP=On to the cmake build command.")
#  endif()
#  message(STATUS "GMP Include: ${GMP_INCLUDE_DIRS}, GMP Library: ${GMP_LIBRARY_DIRS}")
endif()

# Set include directories and link libraries to your project
target_include_directories(MtKaHyPar-Include INTERFACE ${GMP_INCLUDE_DIR})
target_link_libraries(MtKaHyPar-Include INTERFACE ${GMP_LIBRARIES})

...

# Add compile flags that enable warnings
...

# Add support for GNU Multiple Precision)
target_compile_options(MtKaHyPar-BuildFlags INTERFACE -lgmpxx -lgmp)

```

### Include GMP
```cpp
#include <gmpxx.h> // GMP C++ interface
```

## Side trip: k
### Skip some parts if k is large
`LARGE_K` preset skips some allocations to spare time (and memory?). 
We want to skip them too, if `k >= 1024`:

`context.cpp, h`:
```cpp
  bool accountForLargeK() const {
    switch (preset_type) {
      case PresetType::large_k:
        return true;
      case PresetType::cluster:
        return k >= 1024;
    }
    return false;
  }
```
&rarr; for now only used in `simple_rebalancer.cpp` to skip moving nodes to empty parts [Note: `simple_rebalancer` isn't used in `cluster` preset...]

All other places [that I've found] are `switch`-es or skipped `MemoryPool`-registration &rArr; are specific to `LARGE_K` preset...)

### Not necessary to enter blocks (k), epsilon (eps) for clustering
- `command_line_options.cpp`:
    - Disable `.required()` for `-k / --blocks`, `-e / --epsilon` parameters
    - after parsing, throw an error if neither `cluster` preset nor `Objective::aon_hypermodularity` is used:
    ```cpp
    // Validate that blocks is specified
    if ( (context.partition.preset_type != PresetType::cluster &&
          context.partition.objective != Objective::aon_hypermodularity) 
        && !cmd_vm.count("blocks")) {
      throw po::error("The --blocks option is required when the preset is not 'cluster' and objective is not 'aon_hypermodularity'");
    } else {
        context.partition.k = 2;
    }
    // Validate that Epsilon is specified
    if ( (context.partition.preset_type != PresetType::cluster &&
          context.partition.objective != Objective::aon_hypermodularity) 
        && !cmd_vm.count("epsilon")) {
      throw po::error("The --epsilon option is required when the preset is not 'cluster' and objective is not 'aon_hypermodularity'");
    } else {
        context.partition.epsilon = 10000.0; // effectively disables imbalance checks
    }
    ```

## V-Cycles
Per default, IP is skipped in V-cycles and instead the partition from the previous cycle is used as initial partition for the next cycle.

**We want to run IP in each V-cycle, when looking for clusters &rArr; we resompute coefficients of hypermodularity and partitioner again as if we start from scratch with the previous partition as community detection.**

To save memory, we fit the new communities s.t. they use consecutive IDs from `0` to `new_k - 1`. This is done in `fitCommunityIDs()` in time $O(#hn)$ (worst case of one thread). The technic is the same as used in `StaticHypergraph::contract(..)` to renumber the nodes of the contracted hypergraph.

The initial `k` value is set to the number of detected communities. `c-t = 20` (for now) &rarr; 1, 2, 5 return less clusters.

**!!!** In the library interface, `mt_kahypar_partition(..)` sets `num_vcycles = 0`. To use V-Cycles, `mt_kahypar_improve_partition(..)` must be called with `num_vcycles > 1`.

`mt-kahypar/partition/multilevel.cpp`:
- in `multilevel_partitioning(..)`:
```cpp
    ////////////// The k value can be changed here 
    // Only with clustering, singleton sets k = num_nodes !!!
    PartitionID new_k = context.partition.k; 

    if (is_vcycle && context.partition.clustering) {
      timer.start_timer("preprocessing", "Preprocessing");
      new_k = hypergraph.fitCommunityIDs();
      if (context.usesHypermodularityIP()) {
        hypergraph.computeAONParameters();
      }
      timer.stop_timer("preprocessing");
      std::cout << "V-Cycle: fit k to " << new_k << std::endl;
    }

    //////////////////////////////// Change k (1/4)
    if (new_k != context.partition.k && new_k > 1) {
      context.partition.k = new_k;
      timer.start_timer("preprocessing", "Preprocessing");
      context.setupPartWeights(hypergraph.totalWeight());
      context.setupContractionLimit(hypergraph.totalWeight());
      context.setupThreadsPerFlowSearch();
      timer.stop_timer("preprocessing");
    }
    /////////////////////////// End of changing k (1/4)
    
    ...
    /////////////////////////// End of changing k (2/4)

    if ( !is_vcycle || context.partition.clustering) {
        /// Run Initial Partitioning
    } else {
        /// Use partition from previous cycle
    }
```
- in all `PartitionID UnderlyingHypergraph::fitCommunityIDs()` [with slightly different variable names]:
```cpp
  // ! Fit community ids to be consecutive numbers
  PartitionID fitCommunityIDs() {
    Array<size_t> mapping(_num_hypernodes, 0);
    doParallelForAllNodes([&] (const HypernodeID u) {
      mapping[communityID(u)] = 1;
    });
    parallel::TBBPrefixSum<size_t, Array> mapping_prefix_sum(mapping);
    tbb::parallel_scan(tbb::blocked_range<size_t>(UL(0), _num_hypernodes), mapping_prefix_sum);
    PartitionID k = mapping_prefix_sum.total_sum();

    // Remap community ids
    doParallelForAllNodes([&] (const HypernodeID hn) {
      _community_ids[hn] = mapping_prefix_sum[_community_ids[hn]];
    });
    return k;
  }
```

## Clustering
Instead of checking if the clustering preset is used via `context.partition.preset_type == PresetType::cluster` (in `partitioner.cpp`, `multilevel.cpp`, `metrics.cpp`), I introduced `context.partition.clustering` boolean variable.

This way, we can make new presets that use clustering.

## Scripts for experiments
Folder: `_experimental_results/`
- `survey_benchmark/` - all benchmarks from the survey paper [Comparison of modularity-based approaches for nodes clustering in hypergraphs](https://arxiv.org/pdf/2401.14028)


- `build_exp.sh`, `run_exp.sh`, `analyze_exp.sh` - scripts to build, run and analyze the experiments
- `experiments.sh` - uses all 3 scripts above to run many experiments
- `run_full_experiments.sh` - runs `experiments.sh` for many presets, modes and scenarios, plots the results. Uses one of prompts from `prompts/` to interact with `experiments.sh`

### Build
To build and use the library:
- build the library after building the project via: `make install-mtkahypar` after normal `cmake` call
- the `.so` file should be in `LD_LIBRARY_PATH` (e.g. `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/mt-kahypar/lib`)
- in `.cc`: `#include "<mtkahypar.h>` 
- link:  \
    ```g++ -DNDEBUG -o exp -O3 -std=c++20 run_experiments_write_clusters_times_conductances.cc $LIB_DIR/libmtkahypar.so -I .../IFP/include -lpthread```

### Running experiments:
- `run_experiments_write_clusters_times_conductances.cc` [partiton.out] - runs the experiments and writes the resulting clusters, running times and conductances:
    - command line arguments:
        + `-p` - preset (e.g. `dummy` for the preset `dummy_preset.ini` in `config/`)
        + `-v` - version, to save the result in a version subdirectory (`test / ...`)
        + `-m` - name of the mode (`DCHSBM / HyperSBM / HyperSBM`)
        + `-s` - scenario name (`scenA1 / ...`)
        + `-t` - number of threads for the experiment
        + `-instances` - optional: comma-separated list of instance ids to run (default: all instances in the scenario: `1,...,25`)
    - output files:
        + `survey_benchmark/<mode>/<scenario>/<version>/<hgrname>.part` - resulting clusters
        + `survey_benchmark/<mode>/<scenario>/<version>.results.t_c` - running times and conductances:\
            ```
                Version = <version>
                Instances = [i_1, i_2, ..., i_n]          # instance ids
                Time (sec.) = [t1, t2, ..., tn]
                Conductance = [c1, c2, ..., cn]
            ```
- `AON_HMLL_exp.jl` - analog to `analyze_results_via_modularity.jl` (see later)

### Analyzing results:
- `analyze_results_via_modularity.jl` - computes AON-Modularity of the results (and generates results in the case of `AON_HMLL_exp.jl`) and the ground truth clusters and writes a summary file:
    - reference: `AON_HMLL_sript.jl` from [the survey repository](https://github.com/veronicapoda/modularity/)
    - command line arguments:
        + `-d` - data root directory (`survey_benchmark/`)
        + `-m` - name of the generating model (`"DCHSBM/" / "HyperSBM/" / "HyperSBM/"`)
        + `-s` - scenario name (`scenA1 / ...`)
        + `-n` - optional number of instances to analyze (default: all instances in the scenario: 25)
        + `-t` - number of threads for the experiment
        + `-i` - optional: comma-separated list of instance ids to analyze (default: all instances in the scenario)
        + `-v` - version, to read the results from a version subdirectory (`test / ...`)
    - output files:
        + `survey_benchmark/<mode>/<scenario>/<version>.results` (and `<..>.partial_results`) - summary of all stats + previously saved of running times and conductances: \
            ```
                CPU time = [t_1, t_2, ..., t_n]
                Modularity = [m_1, m_2, ..., m_n]
                GT_Modularity = [gtm_1, gtm_2, ..., gtm_m]       
                Relative Modularity Error = [rm_1, rm_2, ..., rm_n] # (GT_Modularity - Modularity) / |GT_Modularity|
                K = [k'_1, k'_2, ..., k'_n]
                GT_K = [k_1, k_2, ..., k_m]             # ground truth number of clusters
                Relative k Error = [rk_1, rk_2, ..., rk_n]       # (GT_K - K) / GT_K
                Instances = [i_1, i_2, ..., i_n]          # instance ids
            ```
- `analyze_results_via_metrics.py` - computes various metrics [ARI, NMI, Purity, F1 Score (Pairwise)] of the ground truth clusters and the given partition and writes a summary file:
    - command line arguments:
        + `--version` - version, to save the result in a version subdirectory (`test / ...`)
        + `--data_root_dir` - data root directory (`survey_benchmark`)
        + `--mode` - name of the generating model (`"DCHSBM/" / "HyperSBM/" / "HyperSBM/"`)
        + `--scenario` - scenario name (`scenA1 / ...`)
        + `--graphs` - optional comma-separated list of instances to analyze
    - output files [APPENDS, if exists!]:
        + `survey_benchmark/<mode>/<scenario>/<version>.results` - summary of all stats : \
            ```
            NMI = [n_1, n_2, ..., n_n]                  # normalized mutual information (see [wikipedia](https://en.wikipedia.org/wiki/Mutual_information#Normalized_variants))
            Purity = [p_1, p_2, ..., p_n]               # purity (see [wikipedia](https://en.wikipedia.org/wiki/Purity_(clustering)))
            ...
            ```
- `analyze_results_via_conductance.cc` [`conductance.out`]- computes conductance of the ground truth clusters and the given partition and writes a summary file:
    - command line arguments:
        + `-p` - preset (e.g. `dummy` for the preset `dummy_preset.ini` in `config/`)
        + `-v` - version, to save the result in a version subdirectory (`test / ...`)
        + `-m` - name of the generating model (`"DCHSBM/" / "HyperSBM/" / "HyperSBM/"`)
        + `-s` - scenario name (`scenA1 / ...`)
        + `-instances` - optional comma-separated list of instances to analyze (default: all instances in the scenario)
    - output files [APPENDS!]:
        + `survey_benchmark/<mode>/<scenario>/<version>.results` - summary of all stats : \
            ```
            ...           
            Conductance = [c'_1, c'_2, ..., c'_n]
            GT_Conductance = [c1, c2, ..., cn]       # ground truth conductances
            ```
### Create new experiment presets
- `_experimental_results/config/` - folder with experiment presets:
    - `dummy_preset.ini` - a dummy preset to create new from
    - `generate_presets.sh` - script to generate new presets from `dummy_preset.ini`

### Plotting results
- `plot_all_results.py` - plots all results from given versions (e.g. `test,AON_HMLL_clique`):
    - command line arguments:
        + `--data_root_dir` - data root directory (`survey_benchmark`)
        + `--versions` - comma-separated list of versions, which results should be plotted (e.g. `test,AON_HMLL_clique`)
        + `--output` - directory to save the plots (e.g. with `--output=AON` the script will generate plots in `plots/AON/`)
