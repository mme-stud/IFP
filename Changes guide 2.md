# Changes Guide: SoSe25

## Ideas:
- in `computeAONParameters`, use `nodeWeightedDegree` for `ClosVol` instead of `nodeDegree` \
(as cutting edges are considered vith weights) and `totalVolume` instead of `initialTotalVertexDegree` for `vol_H` &rarr; done (for now)

## TODO:
- TODO: ensure, that `use_community_detection` is enabled by `aon_hypermodularity` IP [`partitioner.cpp preprocess(..)`] 

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
    Explicitly, $\_beta[k] = - \beta_k$, $\_gamma[k] = - \beta_k * \gamma_k$
    ```cpp
        // AON HyperModularity Clustering Coefficients
        vec<double> _beta;                 ///< -β_k
        vec<double> _gamma;                ///< -β_k * γ_k
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
- getters in `partitioned_hypergraph.h` and a mirroring interface in `partitioned_graph.h`


### Implement AON-Hypermodularity IP

- `aon_hypermodularity_initial_partitioner.h, .cpp`:
    - "Template": `mt-kahypar/partition/initial_partitioning/random_initial_partitioner.h, .cpp`
Plan:
0. [done] maintain original edge size in Static Hypergraph + **mirroring public interface in partitioned** / static / dynamic (hyper)graphs:
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
            void snapshotOriginalEdgeSizes(); // saves sizes
            HypernodeID originalEdgeSize(HyperedgeID e) const;
            HypernodeID originalMaxEdgeSize(HyperedgeID e) const; 
        private:
            HypernodeID _original_max_edge_size; // set at the moment of snapshot
    };
```
1. coarsest_underlying_hg should forget initial weighted degrees etc (as otherwise the original edge sizes are too bog => delta to slow)
   -> use factory to copy the hypergraph (check before doing!)
2. contract static hg, make partitioned from it = collapse
3. after collapse save mapping to map: 
	for new_hn: map[community_id[new_hn]] = new_hn
   before each move update map:
	map[community_id[hn]] = map[community_id[hn_of_new_label]]
	(hn_of_new_label should be availible via he in A_i)
4. at expand adjust z:
	for hn in H:
		z[hn] = map[z[hn]]

### Introduce of the new IP to the framework

### Problems