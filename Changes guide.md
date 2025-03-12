# Guide to my code changes
### Potential problems:
- uncontraction in relative contraction order?: reversed 
    &rarr; should be ok
- implemented only for dynamic hypergraph: for other (graphs... &rarr; `lib_generic_impls.h`)? 
    &rarr; **TODO** for static
- `removeIncidentNets(...)`: changes head of IncidentNetArray &rArr; not thread safe? \
    (was so before &rarr; ok)
- `enable(he)`, `disable(he)` - only in remove single / parallel? Otherwise should weighted deg be updated? \
	(deg is not updated &rarr; should be ok)
- meaning of node weights?  &rarr; arbitrary, default: 1 at start &rArr; no implications should be done
- `setNodePart`, `changeNodePart`: should be thread safe? 
	&rarr; TODO later (?) locks for edges on change pin_count/etc ?..
- `freeInternalData()` [`partitioned_hg.h`]: should I clear `part_volumes`, `part_cut_weights?` \
	(part_weights not cleared &rArr; no) 
- what is `tbb::static_partitioner()`? Is it needed for `initializePartCutWeights()` in `partitioned_hd`?
	&rarr; mode of part algo that we use

### TODO: 
- change `sync_update` to contain info for calculating concuctance gain (1-st and 2-nd max cond cuts; total val; ...)
    - &rarr; **Problem**: conductance gain depends on both whole parts, not just the edge

### Next steps:
- implement `_total_volume`, weighted degrees + "their" functions in `static_hypergraph.h/.cpp`
- write tests (better somewhere in-between)
- implement pq for all conductances (pairs `<Vi_vol, Vi_cut_weight>`)
	- **!!!** Problem: conductance depends on `total_vol` \
	&rarr; heap will not be correct after changes in `total_val` :( \
	    **???** (= contract, uncontract, remove single-pin?) \
	&rarr; Naive: rebuild heap after changes in `total_val`? \
	**???** &rarr; implement `increase_total_val` and refactor to do smth on changes on `total_val`..?
- follow the guide to implement a custom objective function \
	&rarr; adjust `sync_update` (see TODO above)
- debug (potentially many times :( )
...
- outlook: implement new Partitioning Configurations (`include/lib_helper_functions.h`)

## Part 1: needed information
For hypergraph:
+ total volume
+ weighted degrees of \
&rarr; done for `dynamic_hypergraph`, needed for `static_hypergraph`

For partitioned hypergraph:
+ access to total volume, weighted degrees
+ number of pins of e in Vi
+ volume of Vi
+ cut weight of Vi

### Part 1.1: dynamic hypergraph (total volume, weighted degrees)
**???** Also implement for other (heper-)graphs (not only dynamic hg?) as in `lib_generic_impls.h`?
&rarr; for static hg

**Total Volume**:
`dynamic_hypergraph.h`:
+ \+ `total_volume` - `std::atomic<HyperedgeWeight>` due to changes (analog. to `_contraction_index`)
+ [analog. to `_total_weight, totalWeight`]
+ `removeEdge(he)`, `removeLargeEdge(he)`: update `_total_volume`
+ `restoreLargeEdge(he)`: update `_total_volume` \
**???** Why no update for `_total_weight` on `enable()`, `disable()`, `setNodeWeight()` ... ? \
**!!!** &rarr; I do no updates for `_total_weight` (analog.)

`dynamic_hypergraph.cpp`:
- \+ updateTotalVolume()
- analog. to `updateTotalWeight()` !!! but edges instead of nodes !!	
- \+ copy: `_total_volume` (analog. to `_total_weight`)
- `restoreSinglePinAndParallelNets(..)`: update `_total_volume` for single-pin nets \
    **!!!** parallel &rArr; weights are added \
    **!!!**	&rArr; no update of `_total_volume` (and `weighted_degrees`) needed
- `removeSinglePinAndParallelHyperedges()`: update `_total_volume`
- `uncontract(..)`: update `_total_volume` \
    **!!!**&rarr; Important: (non)shared he \
	    (`weighted_degree` is handled by `incident_net_array` &rarr; ok) \
	    TODO (what? &rarr; look there) 
- `contract(v, ..)` &rarr; `contract(u, v, ..)`
    **!!!**&rarr; in `contractHyperedge(u, v, he)`: update `_total_volume`
			(shared he's are removed from `incident_nets` of `v`)

&rarr; adjust `dynamic_hypergraph_factory.cpp`:
- `construct(...)`: compute `_total_volume`

&rarr; adjust `dynamic_hypergraph_test.cc`:
    [hypergraph: from `hypergraph_fixture.h`]
- `HasCorrectStats`, `ModifiesEdgeWeight`, `ComparesStatsIfCopiedParallel`, `ComparesStatsIfCopiedSequential`


**Weighted Degrees (!!!)**
(**???** `dynamic_adjacency_array.cpp` - for vectors **???**)

&rarr; `dynamic_hypergraph.h`:
- \+ `nodeWeightedDegree(u)` analog. to `nodeDegree()` 
    **!!!** `decreaseNodeWeightedDegree(u, w)`: for dealing with single-pin he (`remove`, `restore`)

&rarr; `incident_net_array.h'
- **!!!** \+ `#include .../dynamic_hypergraph.h` 
- \+ `_hypergraph_ptr` (to compute weighted degrees, **!!!** no support for `set/changeWeight`, can be `nullptr`)
    - &rarr; changes in 2 constructors + usages of constructors to include `_hypergraph_ptr` 
	- `dynamic_hypergraph_factory.cpp`: in `construct(...)`
	- `incident_net_array.cpp`: in 2 `copy(...)`-s
- \+ `nodeWeightedDegree(u)` 
- **!!!**  \+ `decreaseNodeWeightedDegree(u, w)`: for dealing with single-pin he (`remove...`, `restore...`) \
	in `Header` analog. to `IncidentNetArray::Header::degree`

&rarr; `incident_net_array.cpp`
- `contract(u, v, ...)`: add `weighted_degree`
- `uncontract(u, v, ...)`: subtract `weighted_degree`  
- `removeIncidentNets(u, ...)`: subtract `weighted_degree` 
    - **!!!** \+ `bool update_weighted_degrees = true` option for correct `removeSinglePinAndParallelHyperedges`
- `restoreIncidentNets(u, ...)`: add `weighted_degree`
    - **!!!** \+ `bool update_weighted_degrees = true` option for correct `removeSinglePinAndParallelHyperedges`
- `construct(...)`: calculate `weighted_degree` analog. to degree (`.local()` etc) \
    **!!!** for degree,  `ThreadLocalCounter = tbb::enumerable_thread_specific< parallel::scalable_vector< size_t >` is used \
    &rarr; for weighted degree, I use `tbb::enumerable_thread_specific< parallel::scalable_vector<HyperedgeWeight> >` \
    [analog. to `IncidentNetArray::Header::degree`]


### Part 1.2: partitioned hypergraph (access to hgInfo, number of pins, volumes, cut weights)
Access to new hypergraph infos:
&rarr; `partitioned_hypergraph.h`:
- \+ `totalVolume()` analog. to `totalWeight()`
- \+ `nodeWeighteddegree(u)` analog. to `nodeDegree(u)`

**Number of pins**: 
- &rarr; `partitioned_hypergpaph.h`: `pinCountInPart(he, id)` &rarr; OK

**Volumes:**
&rarr; `partitioned_hypergpaph.h`
- \+ `_part_volumes`: (vec<CAtomic> **!!!**)
- constructors, `resetData()` - trivially adjusted analog. to `_part_weights`
- \+ `decrementVolumeOfBlock(p, w)`, `incrementVolumeOfBlock(p, w)` \
        analog. to `incrementPinCountOfBlock(e, p)`
- \+ `recomputePartVolumes()` analog. to `recomputePartWeights()`
- \+ `applyPartVolumeUpdates(p_v_deltas)`: for `initializeBlockVolumes` \
        analog. to `applyPartWeightUpdates(p_w_deltas)`
- \+ `initializeBlockVolumes()` - analog. `initializeBlockWeights()` \
    **!!!** uses `nodeWeightedDegree()` &rarr; smth might be wrong with `enable(he)` :(
- &rarr; `initializePartition()`: calls `initializeBlockVolumes()` parallel to weight, pin count
- **!!!** `setEdgeWeight(e, w)`, `enableHyperedge`, `disableHyperedge` - not touched for now
- `uncontract(batch, ...)` update `_part_volumes` (analog. to `_total_volume`)
- `restoreLargeEdge(he)`: update `_part_volumes` analog. to `ets_pin_count_in_part` (`.local()`...)
- `restoreSinglePinAndParallelNets(..)`: update `_part_volumes` for each single-pin net
    - **!!!**&rarr; no need for "thead-safe" behavior: atomic `_part_volumes`
- `setOnlyNodePart(u, p)`: not touched (sets **only** partId)
- `setNodePart(u, p)`: [used *only* to set the first partId of a node] \
	update `_part_volume[p] += nodeWeightedDegree(u)`
- `changeNodePart(u, from, to, ...)`: update `_part_volumes` for `from` and `to`
- **???** `freeInternalData()` : used by destructor for external memory usages (connectivity...) \
	&rArr; nothing changed


**Cut weight of Vi**:
&rarr; `partitioned_hypergpaph.h`
- \+ `_part_cut_weights` (vec<CAtomic> **!!!**)
- constructors, `resetData()` - trivially adjusted analog. to `_part_weights`
- \+ `decrementCutWeightOfBlock(p, w)`, `incrementCutWeightOfBlock(p, w)` \
	analog. to `incrementPinCountOfBlock(e, p)` \
	**???**	potentially reimplement later in ~ `pin_count.h`?
- \+ `recomputePartCutWeights()` analog. to `recomputePartWeights()` ('for testing') 
- \+ `applyPartCutWeightUpdates(p_c_w_deltas)`: for `initializeBlockCutWeights`
	analog. to `applyPartWeightUpdates(p_w_deltas)`
- \+ `initializeBlockCutWeights()` - analog. to `initializeBlockWeights()` \
	**!!!** uses `connectivity(he)` &rarr; run only after `initializePinCountInPart()` \
	**???**	what is `tbb::static_partitioner()?` &rarr; part. algo. that uses `StaticHypergraph`
- &rarr; `initializePartition()` calls `initializeBlockCutWeights()` after all other 
- **!!!** `setEdgeWeight(e, w)`, `enableHyperedge`, `disableHyperedge` - not touched for now
- `restoreLargeEdge(he)`: update `_part_cut_weights`, if `he` - cutting edge
- `setOnlyNodePart(u, p)`: not touched (sets *only* partId)
- `setNodePart(u, p)`: [used only to set the first partId of a node]
	- update `_part_cut_weights` for incident he-s \
		&rarr; if he - new cutting edge (connectivity: 1 &rarr; 2): \
			for both 2 part
		&rarr; if p - new part for cutting edge he (con.: x &rarr; x + 1 | x != 1): \
			only for part p 
	- **!!!** &rArr; now uses `connectivity_info`, `connectivity_set` \
	&rArr; should be run **sequentialy** 
- `changeNodePart(u, from, to, ...)` uses `updatePinCountOfHyperedge(he,from,to,sync_update,..)`
- &rarr; `updatePinCountOfHyperedge(he, from, to, sync_update, ..)`: \
	uses `_pin_count_update_ownership[he]` lock to ensure thread-safety! \
	&rarr; update `_part_cut_weights` for `from` and `to`
- **???** `freeInternalData()` : used by destructor for external memory usages (connectivity...) \
	&rArr; nothing changed