# Changes Guide: SoSe25

## Other changes, Ideas:
- in `computeAONParameters`, use `nodeWeightedDegree` for `ClosVol` instead of `nodeDegree` \
(as cutting edges are considered vith weights) and `totalVolume` instead of `initialTotalVertexDegree` for `vol_H` &rarr; done (for now)
- ~~in `partitioner.cpp` change `k` to the number of found communities, if `cluster && community_detection=true` (for now, I just recalculate the number of communities as done in `partitioner_output.cpp`)~~ [this way the kernel is bigger &rArr; AON IP couldn't find a good clustering and returned too many clusters...]
- in `multilevel.cpp` before uncoarsening stage check if k has changed (+correctly change it, if so) &rarr; done, not used for now
- **!!! I concider edge weights in the gain &rArr; use weighted degrees and use edge weight in _delta_cut** &rarr; ASK
- enable single pin net removal (they are never cutting and volumes are tracked correctly)
- [done] set `_beta[d] = 0` if `NaN`
- [done] disable balancing for clustering in lp
- set lp to sequential (as it is bad in parallel :( )
- run both original sized and reduced sized versions of AON IP + a couple of random IP's


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
- merge with Adil to allow `double` gains &rarr; maybe than return to `delta` = sum of all local "gains" by `conductance_local` in `refineImpl()` (`label_propagation_refiner.cpp`) 
- try to improve local search (gain cache? another approach?)


## Initial Partitioning: Hypermodularity

*(new!)* `aon_hypermodularity` uses original edge sizes (if possible), `aon_hypermodularity_kernel` uses the edge sizes in the kernel.

ToDo:
- calculate parameters ($\omega_{k0}, \omega_{k1} \beta, \gamma$) needed for the AON-Hypermodularity IP &rarr; community detection in preprocessing
- implement Hypermodularity as an IP (no parallelization in IP as multiple IP algorithms could be ran simultaneously + kermel should normally be small)
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
0. use `long double` in volume calculations to avoid `Inf` as much as possible (withput using GMP...) *[Note: if `long doubles` are used for beta and gamma too, then the running time increases significantly (from <1 min to 10-20 min)]*
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
                     const bool randomize,
                     bool useOriginalEdgeSizes);
  
  void partitionImpl() final {
    partitionImpl(
      1e3 /* edgeSizeThreshold */,
      1e2 /* maxNumIter */,
      1e-8 /* eps */,
      true /* randomize */,
      true /* useOriginalEdgeSizes */
    );
  }
```
- Main code:
```cpp
      if (! useOriginalEdgeSizes || ! H.hasOriginalEdgeSizes()) {
        useOriginalEdgeSizes = false;
        if ( ! H.hasOriginalEdgeSizes() )
          LOG << "AON IP: No snapshot of original edge size found";
        // save current edge sizes, weighted degrees and total volume
        H.snapshotOriginalEdgeSizes();
        H.snapshotOriginalWeightedDegreesAndTotalVolume();
        H.useOriginalSizeInParallelNetsDetection(true); // otherwise gain is incorrect
      }
      if (useOriginalEdgeSizes)
        LOG << "AON IP: Using original edge sizes.";
      else
        LOG << "AON IP: Using current edge sizes.";
      H.enableSinglePinNetsRemoval(); // single pin nets are never cutting

      //          1. Singleton initial partitioning
      //                         <...>
      //          2. AllOrNothingHMLL: Louvain Cycle   
      double total_gain = 0.0;
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
        collapse(H_new, H_new_partitioned, map_z);

        /** ------------------ Louvain Step: ------------------
         *  - Nodes are moved to neighbouring partitions as
         *    long as it improves the modularity gain;
         *  - map_z is updated accordingly.
         */
        new_gain += louvainStep(H_new, H_new_partitioned, map_z, _beta, _gamma, 
                                  edgeSizeThreshold, maxNumIter, eps, randomize);
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
      }

      LOG << "AON IP finished Louvain with total gain " << total_gain;

      //             3. Finalize Partitioning
      //                        <...>

```

`aon_hypermodularity_kernel_initial_partitioner.h, .cpp`:
- a single attribute `AONHypermodularityInitialPartitioner _aon_ip`
- calls `aon_hypermodularity` with `useOriginalEdgeSizes=false`

#### Needed Additional Functionality (Static Hypergraph, Partitioned Hypergraph)

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
      LOG << "PartitionedHypergraph::fitK() - Warning: only one cluster found: "
             "actualK = " << actualK << ", setting it to 2";
    }
    if (_k != actualK) {
      setK(actualK);
      LOG << "PartitionedHypergraph::fitK() - Fitted k to " << actualK;
    }
  }

```



#### Implementation Details
0. The underlying hypergraph `H` can have too many single-pin hyperedges \ 
&rArr; I contract it's singleton communities and after that disabled single-pin nets removal
1. Louvain `Collapse(..)` and `Expand(..)`:
    ![Algorithm 3](<Algorithm 3: AllOrNothingHMLL.png>)
    - **Main idea**: contract static hg (*= the result of the last contraction*) and create a partitioned hg from it = collapse
    - after `collapse(..)` save mapping to `map_z`: \
    	`map_z[communityID(collapsed_hn)] = collapsed_hn`
    - at `expand(..)` adjust `z` (*= the best found partitioning of the given coarsest hypergraph*): \
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
        while (improving && (iter++ < maxNumIter)) {
            improving = false;
            double gain = 0.0;
            if (randomize) {
                vec<HypernodeID> nodes(numNodes, 0);
                for (HypernodeID i = 0; i < numNodes; ++i) {
                    nodes[i] = i;
                }
                std::shuffle(nodes.begin(), nodes.end(), _rng);
                for (const HypernodeID &i : nodes) {
                    gain += louvainStepForANode(i, neighbours[i], visited, H_new_partitioned, map_z, beta, gamma, maxNumIter, eps, randomize, edgeSizeThreshold);
                }
            } else {
                for (const HypernodeID &i : H_new_partitioned.nodes()) {
                    gain += louvainStepForANode(i, neighbours[i], visited, H_new_partitioned, map_z, beta, gamma, maxNumIter, eps, randomize, edgeSizeThreshold);
                }
            }
            total_gain += gain;
            improving = (gain > eps);
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

    [For now] I introduced `const HypernodeID edgeSizeThreshold=1e3` parameter to avoid too big loops in gain computation and (more importantly) calculating with infinities (as inf - Inf = NaN &rarr; the gain is `NaN` &rArr; move isn't done)

### Introduce of the new IP to the framework

1. Adjust the parameters of `cluster` preset:
    - set `k` to 32 instead of 2 (Adil has done so -> TODO: ask):
        - `mt-kahypar\partition\context.cpp` in `Context::setupContractionLimit(total_hypergraph_weight)`: \
        for `cluster` preset, set `coarsening.contraction_limit` to `coarsening.contraction_limit_multiplier * 32` instead of `coarsening.contraction_limit_multiplier * partition.k`.
        - `mt-kahypar\partition\partitioner.cpp` in `setupContext(& hypergraph, & context, *target_graph)`: \ 
        set `k=32` for `cluster` preset
    - [old changes guide] analog. to `singleton` introduce `aon_hypermodularity`, `aon_hypermodularity_kernel`:
        - in `config/`:
            - `cluster_preset.ini`: 
                - set `i-enabled-ip-algos=1` only for `aon_hypermodularity` and `aon_hypermodularity_kernel`
            - `large_k_preset.ini`: 
                \+ `i-enabled-ip-algos=0` for `aon_hypermodularity`, `aon_hypermodularity_kernel`
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
                    UNDEFINED = 12
                    };
                    ```
                - `context_enum_classes.cpp`:
                    - in `operator<<(or, algo)` and `initialPartitioningAlgorithmFromString(algo)` add transmations string <-> `InitialPartitioningAlgorithm` for `aon_hypermodularity`, `aon_hypermodularity_kernel`
                - in `initial_partitioning/`: 
                    - \+ `aon_hypermodularity_initial_partitioner.h, .cpp`
                    - \+ `aon_hypermodularity_kernel_initial_partitioner.h, .cpp` 
                    - **!!!** `CMakeLists.txt`: \+ `aon_hypermodularity_initial_partitioner.cpp`, `aon_hypermodularity_kernel_initial_partitioner.cpp`
                - in `registries/`:
                    - `register_initial_partitioning_algorithms.cpp`
                    - \+ `#include "../initial_partitioning/aon_hypermodularity_initial_partitioner.h"`, `#include "../initial_partitioning/aon_hypermodularity_kernel_initial_partitioner.h"`
                    - \+ define `AONHypermodularityPartitionerDispatcher`, `AONHypermodularityKernelPartitionerDispatcher`
                    - in `register_initial_partitioning_algorithms()`: \+ register `AONHypermodularityPartitionerDispatcher`, `AONHypermodularityKernelPartitionerDispatcher`
            - in `io/`:
                - `command_line_options.cpp`: 
                    - by `"i-enabled-ip-algos"` example add aon_hypermodularity, aon_hypermodularity_kernel IP (and change the number of IP-algos at the end of the example)
                - `presets.cpp`:
                    - in `load_large_k_preset()`: by`// main -> initial_partitioning` add entry for `aon_hypermodularity` (`"0"`), entry for `aon_hypermodularity_kernel` (`"0"`)
                    - in `load_clustering_preset()`: 
                        ```cpp
                         "0" // singleton" IP
                         "1" // aon_hypermodularity
                         "0" // aon_hypermodularity_kernel (always worse)
                        ```
    - `config/cluster_preset.ini` and `mt-kahypar/io/presets.cpp`: in `# main -> initial_partitioning` set `i-runs=3` istead of 10 for `cluster` (as AON-hypermodularity is ~~deterministic~~ randomized and runs `i-runs * t` times)

3. `sanity_check(*target_graph)` in `context.cpp`:
    - adjust conductance checks to allow `aon_hypermodularity`, `aon_hypermodularity_kernel` IP
    - ensure, that `use_community_detection` is enabled if a hypermodularity IP is used

4. ~~[Idea] Change `context.partition.k` if it changed after IP (due to `aon_hypernodularity`)  &rarr; not done, as `new_k` shouldn't be greater than the number of nodes~~
4. change `context.partition.k` in `multilevel.cpp` if it `aon_hypernodularity` IP is used \ 
[analog. to `cluster` + `singleton` &rArr; `new_k = #_nodes`] 

#### Define helper functions to detect hypermodularity IP

`mt-kahypar/partition/context_enum_classes.h, cpp`:
- \+ `bool isHypermodularityIP(InitialPartitioningAlgorithm algo)`

`mt-kahypar/partition/context.h`:
- \+ `bool usesHypermodularityIP() const`

### Run Hypermodularity IPs on parallel
Hypermodularity IP heavily depends on the visiting order of the nodes.
&rArr; we run multiple hypermodularity IPs in parallel.

**For now for hypermodularity IPs, I multiply the number of runs by the number of threads.**

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
- disable error of missing `-k / --blocks`, `-e / --epsilon` parameter for `aon_hypermodularity` (and `conductance_local`, `conductance_global`):
    - `mt-kahypar/io/command_line_options.cpp`: `processCommandLineInput(...)`
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

## Label Propagation

### Problems:
- more than 3 threads &rArr; bad results :(
    - **ToDo (?)**: set a sequential flag for LP to be able to run multiple IP in parallel
- [solved for now] lp refiner summs up all Attributed gains of the moves in a `delta`. With `HyperedgeWeight` gains I get an overflow &rArr; `refineImpl(..)` returns `detla > 0  = false` and the refinement round finishes:
    - Current solution: recalculate `delta` with the actual conductance gain
    - To be changed (?) after introducting `double` gains
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
    - **ToDo**: implement a (dummy?) gain cache