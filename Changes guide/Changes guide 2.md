# Changes Guide: SoSe25

## Other changes, Ideas:
- in `computeAONParameters`, use `nodeWeightedDegree` for `ClosVol` instead of `nodeDegree` \
(as cutting edges are considered vith weights) and `totalVolume` instead of `initialTotalVertexDegree` for `vol_H` &rarr; done (for now)
- in `partitioner.cpp` change `k` to the number of found communities, if `cluster && community_detection=true` (for now, I just recalculate the number of communities as done in `partitioner_output.cpp`)
- in `multilevel.cpp` before uncoarsening stage check if k has changed (+correctly change it, if so) &rarr; done, not used for now
- **!!! I concider edge weights in the gain &rArr; use weighted degrees and use edge weight in _delta_cut** &rarr; ASK
- enable single pin net removal (they are never cutting and volumes are tracked correctly)


## Questions
- singleton &rarr; conductance = 1. Why?

## TODO:
- TODO: ensure, that `use_community_detection` is enabled by `aon_hypermodularity` IP [`partitioner.cpp preprocess(..)`] &rarr; done in `context.cpp sanity_check(..)`
- change `context.partition.k` in `multilevel.cpp` if `aon_hypernodularity` IP is used &rArr; done
- introduce ClusteringMode to mark that `k` can be changed (?)
- run multiple AON IP in paralle (#AON = #threads) &rarr; **Problem:** parallel local search works bad: set it to sequential?
- try to snapshot and preserve original edge sizes from the start (if they are not too big: <500?)
- compare AON IP with Hypermodularity (+metrix)
- compare clustering with Hypermodularity (metric = AON &rarr; Adil has implemented the objective)

- iterate between local and global conductance in local search
- merge with Adil to allow `double` gains
- try to improve local search (gain cache? another approach?)


## Initial Partitioning: Hypermodularity

ToDo:
- calculate parameters ($\omega_{k0}, \omega_{k1} \beta, \gamma$) needed for the AON-Hypermodularity IP &rarr; community detection in preprocessing
- implement Hypermodularity as an IP (no parallelization in IP as multiple IP algorithms could be ran simultaneously + kermel should normally be small)
- introduce the new IP to the system:
    - register IP as an IP
    - set it as default IP for `cluster` preset
    - adjust the AON-parameter calculation to take place only if the IP is chosen
- try and debug

### Calculate the parameters
Reference: [commit](https://github.com/adilchhabra/mt-kahypar/commit/ab9be0777bbe77c158bf8e6f53166ea3c67ce526) My change: weighted degrees, total volume.
- `partition\partitioner.cpp`:
    - `precomputeHyperModularityParameters(&hypergraph, &context)` [my: runs only if `aon_hypermodularity` ip algo is enabled] - called by `preprocess(&hg, &context, *target_graph)` if `use_community_detection == true`
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
1. Avoid `NaN` in omegas, by avoiding `vol_out = Inf - Inf = Nan`. If `vol_in = +Inf`, then `vol_out` should be set to `+Inf` as well \ 
    (as $(vol\_H)^d = (\sum_i vol(i))^d, vol\_in = \sum_i vol(i)^d \implies vol\_out = vol\_H^d - vol\_in$, so if `vol_in` is infinite, then `vol_out` should be at least near to `+Inf`)
2. Set `_beta[d] = 0` if it's `Nan` (when both omegas are 0, `_beta[d] = log (+Inf / +Inf) = Inf - Inf = NaN`)
```cpp
    if (std::isfinite(vol_in))
      vol_out = std::pow(vol_H, static_cast<int>(d)) - vol_in;
    else {
      ASSERT(vol_in > 0, "vol_in is not finite, but non-positive: " << vol_in);
      // 1 (Avoiding inf - inf = NaN)
      vol_out = std::numeric_limits<double>::infinity();
    } 

    double omega_in = (m_k[d] - cut_k[d]) / vol_in; 
    double omega_out = cut_k[d] / vol_out;
    // (Normally) not NaN / +-Inf as m_k, cut_k < total_volume. (vol_in and vol_out != NaN, 0)

    _omega[d] = {omega_in, omega_out};
    _beta[d] = std::log(omega_in) -  std::log(omega_out);
    _gamma[d] = omega_in - omega_out;
    
    // Adjustments for Inf, Nan
    if (!std::isfinite(_beta[d])) {
    ASSERT(_beta[d] > 0 || _beta[d] != _beta[d], 
        "_beta[" << d << "] is not finite, not +inf and not NaN: " << _beta[d]);
    _beta[d] = (_beta[d] > 0 ? 1e3 : 0);
    // Idea: _beta[d] = NaN => both log omega-s are the same Inf => _beta[d] = 0
    }
```

#### Problem of huge (original) edges
If the original edge sizes are too big, the `_beta` and `_gamma` vectors are also huge. But *normally* most of their entries are 0.0 (`-1000` and `NaN` before &rarr; see *NaN Problem*).

&rarr; I track the last non-zero pair `(_beta[d], _gamma[d])` and cut all other zeroes. The getters of beta and gamma will return 0.0, if `d` is greater than the last non-zero index.

### Implement AON-Hypermodularity IP
Original Algorithm: [Generative hypergraph clustering: from blockmodels to modularity](https://arxiv.org/pdf/2101.09611)

#### Code Structure 
("Template": `mt-kahypar/partition/initial_partitioning/random_initial_partitioner.h, .cpp`)

- `aon_hypermodularity_initial_partitioner.h, .cpp`:
```cpp
      // save current edge sizes, weighted degrees and total volume
      H.snapshotOriginalEdgeSizes();
      H.snapshotOriginalWeightedDegreesAndTotalVolume();
      H.useOriginalSizeInParallelNetsDetection(true); // otherwise gain is incorrect
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
        total_gain += louvainStep(H_new, H_new_partitioned, map_z, _beta, _gamma, maxNumIter, eps, randomize);

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
      }

      LOG << "AON IP finished Louvain with total gain " << total_gain;

      //             3. Finalize Partitioning
      //                        <...>

```

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
                    gain = louvainStepForANode(i, neighbours[i], visited, H_new_partitioned, map_z, beta, gamma, maxNumIter, eps, randomize, edgeSizeThreshold);
                }
            } else {
                for (const HypernodeID &i : H_new_partitioned.nodes()) {
                    gain = louvainStepForANode(i, neighbours[i], visited, H_new_partitioned, map_z, beta, gamma, maxNumIter, eps, randomize, edgeSizeThreshold);
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
    - [old changes guide] analog. to `singleton` introduce `aon_hypermodularity`:
        - in `config/`:
            - `cluster_preset.ini`: 
                - set `i-enabled-ip-algos=1` only for `aon_hypermodularity`
            - `large_k_preset.ini`: 
                \+ `i-enabled-ip-algos=0` for aon-hypermodularity \
            - (all other `.ini` use `initial_partitioning: i-mode=rb` [recursive bipartitioning] &rArr; no changes)
        - `mt-kahypar/`:
            - in `partition/`:
                - `context_enum_classes.h`:
                    - in `InitialPartitioningAlgorithm`:
                    ```cpp
                    enum class InitialPartitioningAlgorithm : uint8_t {
                    ...
                    aon_hypermodularity = 10,
                    UNDEFINED = 11
                    };
                    ```
                - `context_enum_classes.cpp`:
                    - in `operator<<(or, algo)` and `initialPartitioningAlgorithmFromString(algo)` add transmations string <-> `InitialPartitioningAlgorithm` for `aon_hypermodularity`
                - in `initial_partitioning/`: 
                    - \+ `aon_hypermodularity_initial_partitioner.h, .cpp` :)
                    - **!!!** `CMakeLists.txt`: \+ `aon_hypermodularity_initial_partitioner.cpp`
                - in `registries/`:
                    - `register_initial_partitioning_algorithms.cpp`
                    - \+ `#include "../initial_partitioning/aon_hypermodularity_initial_partitioner.h"`
                    - \+ define `AONHypermodularityPartitionerDispatcher`
                    - in `register_initial_partitioning_algorithms()`: \+ register `AONHypermodularityPartitionerDispatcher`
            - in `io/`:
                - `command_line_options.cpp`: 
                    - by `"i-enabled-ip-algos"` example add aon_hypermodularity IP (and change the number of IP-algos at the end of the example)
                - `presets.cpp`:
                    - in `load_large_k_preset()`: by`// main -> initial_partitioning` add entry for `aon_hypermodularity` (`"0"`) 
                    - in `load_clustering_preset()`: 
                        - `"0" // singleton" IP`, `"1" // aon_hypermodularity`
    - `config/cluster_preset.ini` and `mt-kahypar/io/presets.cpp`: in `# main -> initial_partitioning` set `i-runs=1` istead of 10 for `cluster` (as AON-hypermodularity is deterministic)

3. `sanity_check(*target_graph)` in `context.cpp`:
    - adjust conductance checks to allow `aon_hypermodularity` IP
    - ensure, that `use_community_detection` is enabled if `aon_hypermodularity` IP is used

4. ~~[Idea] Change `context.partition.k` if it changed after IP (due to `aon_hypernodularity`)  &rarr; not done, as `new_k` shouldn't be greater than the number of nodes~~
4. change `context.partition.k` in `multilevel.cpp` if it `aon_hypernodularity` IP is used \ 
[analog. to `cluster` + `singleton` &rArr; `new_k = #_nodes`] 

### Problems