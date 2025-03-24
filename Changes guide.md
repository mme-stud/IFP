# Guide to my code changes
### Notes:
- `HyperedgeWeight` has to be positive (or at least non-negative). Hence, `NonnegativeFraction`.

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
	&rarr; heap will not be correct after changes in `total_vol` :( \
	    **???** (= contract, uncontract, remove single-pin?) \
	&rarr; Naive: rebuild heap after changes in `total_vol`? \
	**???** &rarr; implement `increase_total_vol` and refactor to do smth on changes on `total_vol`..?
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

#### Total Volume:

**For ``DynamicHypergraph``**:

`dynamic_hypergraph.h`:
+ \+ `total_volume` - `std::atomic<HyperedgeWeight>` due to changes (analog. to `_contraction_index`)
+ \+ `totalVolume()`, `updateTotalVolume()`, `updateTotalVolume(par_tag)` : analog. to `_total_weight, totalWeight`
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


**For ``StaticHypergraph``**:
`static_hypergraph.h`:
- \+ `_total_volume`: `HyperedgeWeight` &rarr; think: better unsigned long long?
- \+ `totalVolume()`: analog. to `_total_weight`, `totalWeight`
- `removeEdge(he)`, `remmoveLargeEdge(he)`: update `_total_volume` (used only for testing)
- `restoreLargeEdge(he)`: update `_total_volume` (used only for testing)
- \+ declare `computeAndSetTotalVolume(parallel_tag_t)` analog. to `computeAndSetTotalNodeWeight(parallel_tag_t)`

`static_hypergraph.cpp`:
- `contract(..)`: In Stage 4 [construction of the coarsened hg]: in `tbb::parallel_for` "Write hyperedges from temporary buffers to incidence arrays": recalculate `_total_value` (edge sizes and weights are already correct) \
	**!!!** use thread-local counters as `_total_volue` isn't atomic [debug] - analog. to `_max_edge_size`
- `copy(parallel_tag_t)`, `copy()`: copy `_total_value`
- \+ `computeAndSetTotalVolume(parallel_tag_t)`: analog to `computeAndSetTotalNodeWeight(parallel_tag_t)`. Used by `static_hypergraph_factory.cpp` and `static_hypergraph_test.cc`. **!!!** uses `_weighted_degrees`

&rarr; `static_hypergraph_factory.cpp`:
- `construct(..)`: 
	- compute `_total_volume`: `computeAndSetTotalVolume(parallel_tag_t)` *(! thread-safe !)*

#### Weighted Degrees (!!!)
(**???** `dynamic_adjacency_array.cpp` - for vectors **???**)

**For `DynamicHypergraph`**:

&rarr; `dynamic_hypergraph.h`:
- \+ `nodeWeightedDegree(u)` analog. to `nodeDegree()` 
- \+ **!!!** `decreaseNodeWeightedDegree(u, w)`: for dealing with single-pin he (`remove`, `restore`)
- move constructor, move assigment: call `adjustHypergraphPtr(this)` for `_incident_nets` [debug]

&rarr; `incident_net_array.h'
- \+ `CAtomic<HyperedgeWeight> weighted_degree` of a node in `Header` of its incident_net's list (analog. to `IncidentNetArray::Header::degree`)
	- **!!!**  \+ `decreaseNodeWeightedDegree(u, w)`: for dealing with single-pin he (`remove...`, `restore...`)
	- `CAtomic`, as it is potentially changed simultaniously by `DynamicHypergraph::removeSinglePinAndParallelHyperedges()` through `decreaseNodeWeightedDegree(..)`
- \+ `nodeWeightedDegree(u)` : asserts `_hypergraph_ptr != nullptr`
- **!!!** ~~\+ `#include .../dynamic_hypergraph.h` ~~ forward declaration of `DynamicHypergraph`
- \+ `_hypergraph_ptr` (to compute weighted degrees, **!!!** no support for `set/changeWeight`, can be `nullptr`)
    - &rarr; changes in 2 constructors of `IncidentNetArray` + usages of constructor (?) to pass on `_hypergraph_ptr`:
		- in non-trivial constructor: \+ parameter `const HyperedgeWeight* hyperedge_weight_ptr = nullptr` to pass it in `construct(..)`, to make parallel construction of hg and its `incident_nets` possible (see `construct(..)` in `dynamic_hypergraph_factory.cpp`)
		- \+ new parameter in `construct(..)`: `const HyperedgeWeight* hyperedge_weight_ptr = nullptr` **TODO** check usage
		- **!!!** \+ `public adjustHypergraphPtr(hg_ptr)` to adjust `_hypergraph_ptr` in move constructor, move assigment of `DynamicHypergraph`, (maybe in `DynamicHypergraphFactory` too? &rarr; TODO?)
	- `dynamic_hypergraph_factory.cpp`: in `construct(...)`: \
		pass `const HyperedgeWeight* _hyperedge_weight_ptr` - [const pointer] and `&hypergraph` in `IncidentNetArray(...)` constructor
	- `incident_net_array.h`: in 2 `copy(...)`-s add new parameter `hypergraph_ptr = nullptr`
	- `incident_net_array.cpp`: in 2 `copy(...)`-s ~~copy `_hypergraph_ptr`~~:
		- set `_hypergraph_ptr` of the copy to the passed on `hypergraph_ptr` (passed by 2 `DynamicHypergraph::copy(..)`)
	- &rarr; `dynamic_hypergraph.cpp`:
		- in `copy()` and `copy(parallel_tag)` call `_incident_nets.copy` with a new parameter `&hypergraph` (to set the pointer to the copied hg [**debugging**])

&rarr; `incident_net_array.cpp`
- `#include '.../dynamic_hypergraph.h'` to avoid usage of an incomplete type `DynamicHypergraph` (forward declaration in `incident_net_array.h`)
- `contract(u, v, ...)`: add `weighted_degree`
- `uncontract(u, v, ...)`: subtract `weighted_degree`  
- `removeIncidentNets(u, ...)`: subtract `weighted_degree` 
    - **!!!** \+ `bool update_weighted_degrees = true` option for correct `removeSinglePinAndParallelHyperedges`
- `restoreIncidentNets(u, ...)`: add `weighted_degree`
    - **!!!** \+ `bool update_weighted_degrees = true` option for correct `removeSinglePinAndParallelHyperedges`
- `construct(...)`: calculate `weighted_degree` analog. to degree (`.local()` etc) \
    **!!!** for degree,  `ThreadLocalCounter = tbb::enumerable_thread_specific< parallel::scalable_vector< size_t >` is used \
    &rarr; for weighted degree, I use `tbb::enumerable_thread_specific< parallel::scalable_vector<HyperedgeWeight> >` \
    [analog. to `IncidentNetArray::Header::degree`] \
	**!!!** if `hyperedge_weight_ptr` is passed on from `construct(..)` in `DynamicHypergraphFactory`, we should use it, else use weight=1 for all he instead of `_hypergraph_ptr`, as **the hypergraph is constructed in parallel to its incident net array.** \
	**!!!** Initialize weighted degree to 0 in header(p) as **no `Header` constructor is called** due to pointer tricks with `static_cast` *[debug]*

**For ``StaticHypergraph``**

`static_hypergraph.h`:
- \+ `_weighted_degrees = Array<HyperedgeWeight>`
- \+ `nodeWeightedDegree(u)` trivial 
- \+ `decreaseNodeWeightedDegree(u, w)`: trivial (to ensure same interface)
- `removeEdge(he)`, `remmoveLargeEdge(he)`: update `_weighted_degrees` (used only for testing)
- `restoreLargeEdge(he)`: update `_weighted_degrees` (used only for testing)

`static_hypergraph.cpp`:
- `contract(..)`: Stage 4 [construction of coarsened hg]: resize `_weighted_degrees`, calculate weighted degrees in `setup_hypernodes` \
	&rarr; **LOOK IN `dynamic_hypergraph`, `partitioned_hypergraph`: resized all arrays?** &rarr; yes \
	**!!!** `he.weight()` - only for **enabled** edges &rArr; used `tmp_hyperedges[id].weight()`
- `copy(parallel_tag_t)`: copy `_weighted_degrees` parallel and analog. to other arrays (under `_incidence_array`).
- `copy()`: copy `_weighted_degrees` analog. to other arrays (under `_incidence_array`).
- `memoryConsumption(parent)`: add memory consumption of `_weighted_degrees` array
- STOPPED HERE

&rarr; `static_hypergraph_factory.cpp`:
- `construct(..)`: 
	- resize `_weighted_degrees`
	- compute `_weighted_degrees`: `hyperedge_weight`is a ptr \
		&rArr; `hyperedge_weight ? hyperedge_weight[pos] : 1;` \
		**!!!** use thread-local storage `tbb::enumerable_thread_specific< parallel::scalable_vector < HyperedgeWeight> > local_weighted_degree_per_vertex(num_hypernodes, 0);` analog. to `num_incident_nets_per_vertex` [debug]

### Part 1.2: partitioned hypergraph (access to hgInfo, number of pins, volumes, cut weights)
Access to new hypergraph infos:
&rarr; `partitioned_hypergraph.h`:
- \+ `totalVolume()` analog. to `totalWeight()`
- \+ `nodeWeighteddegree(u)` analog. to `nodeDegree(u)`

#### Number of pins: 
- &rarr; `partitioned_hypergpaph.h`: `pinCountInPart(he, id)` &rarr; OK

#### Volumes:
&rarr; `partitioned_hypergpaph.h`
- \+ `_part_volumes`: (`vec<CAtomic>` **!!!**)
- constructors, `resetData()` - trivially adjusted analog. to `_part_weights`
- \+ `partVolume(p)` - getter
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
- `resetPartition()`: reset `_paer_volumes` analog. to `_part_weights` 
- **???** `freeInternalData()` : used by destructor for external memory usages (connectivity...) \
	&rArr; nothing changed


#### Cut weight of Vi:
&rarr; `partitioned_hypergpaph.h`
- \+ `_part_cut_weights` (`vec<CAtomic>` **!!!**)
- constructors, `resetData()` - trivially adjusted analog. to `_part_weights`
- `partCutWeight(p)` - getter
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
		&rarr; if p - new part for cutting edge he (con.: x &rarr; x + 1 | x + 1 > 1): \
			only for part p 
	- **!!!** &rArr; now uses `connectivity_info`, `connectivity_set` \
	&rArr; should be run **sequentialy** 
- `changeNodePart(u, from, to, ...)` uses `updatePinCountOfHyperedge(he,from,to,sync_update,..)`
- &rarr; `updatePinCountOfHyperedge(he, from, to, sync_update, ..)`: \
	uses `_pin_count_update_ownership[he]` lock to ensure thread-safety! \
	&rarr; update `_part_cut_weights` for `from` and `to`
- `resetPartition()`: reset `_paer_volumes` analog. to `_part_weights`
- **???** `freeInternalData()` : used by destructor for external memory usages (connectivity...) \
	&rArr; nothing changed

### Part 1.3 Test

#### StaticHypergraph, DynamicHypergraph (Total volume, Weighted degrees)
- `dynamic_hypergraph_test.cc`: `ADynamicHypergraph` [check `nodeWeightedDegree(..)`, `totalVolume`]
	- `HasCorrectStats`: `totalVolume == 12`
	- \+ `VerifiesVertexWeightedDegrees`: [analog. to `VerifiesVertexDegrees`] fixture sets no weights (`nullptr`) &rArr; `nodeDegree = nodeWeightedDegree`
	- `ModifiesEdgeWeights`: `updateTotalVolume()`, `totalVolume == 17`
	- `ComparesStatsIfCopiedParallel`, `ComparesStatsIfCopiedSequential`: `totalVolume`
	- `PerformsAContraction1 ... 5`: check both
	- `PerformAContractionsInParallel1 ... 3` check both
	- `verifyEqualityOfDynamicHypergraphs`: check both
	- `RemovesSinglePinAndParallelNets1 ... 2` check both
	- `RestoreSinglePinAndParallelNets1 ... 2` check both
	- `GeneratesACompactifiedHypergraph2` check both
- `hypergraph_fixtures.h`: ~~TODO(?) `verifyPins`: weighted degrees ? &rarr; no!~~ (done in `dynamic_hypergraph_test.cc`)
- `incident_net_array_test.cc`: [check `nodeWeightedDegree(..)`, `_hypergraph_ptr` (?)]
	- TODO tests with _hypergraph_ptr != nullptr (?)
- `static_hypergraph_test.cc`: AStaticHypergraph [check `nodeWeightedDegree(..)`, `totalVolume`]
	- `HasCorrectStats`: `totalVolume() == 12`
	- \+ `VerifiesVertexWeightedDegrees`: `nodeWeightedDegree(..)`
	- `RemovesAHyperedgeFromTheHypergraph1 ... 4`: check both
	- `ComparesStatsIfCopiedParallel`, `ComparesStatsIfCopiedSequential`: `totalVolume`
	- \+ `ComparesWeightedDegreesIfCopiedParallel`, `ComparesWeightedDegreesIfCopiedSequential`: weighted degrees
	- `ContractsCommunities1 ... 3`, `ContractsCommunitiesWithDisabledHypernodes`, `ContractsCommunitiesWithDisabledHyperedges`: check both
 
#### PartitionedHypergraph (PartVolume, PartCutWeight, ConductancePriorityQueue)
- `partitioned_hypergraph_test.cc`: APartitionedHypergraph
	- \+ `HasCorrectPartVolumes`,  analog. to `HasCorrectPartWeightAndSizes` 
	- \+ `HasCorrectPartVolumesIfOnlyOneThreadPerformsModifications`, `HasCorrectPartCutWeightsIfOnlyOneThreadPerformsModifications` analog. to `HasCorrectPartWeightsIfOnlyOneThreadPerformsModifications`
	- `PerformsConcurrentMovesWhereAllSucceed`: check `partVolume`, `partCutWeight`, `conductance_pq`
	- `ComputesPartInfoCorrectIfNodePartsAreSetOnly`: check `partVolume`, `partCutWeight`, `conductance_pq`

- `dynamic_partitioned_hypergraph.cc`: ADynamicPartitionedHypergraph
	- nothing to do (?): only define checking functions for new `Objective`

- `priority_queue_test.cc`: TODO (?) analog. tests for `ConductancePriorityQueue`

- `partitioned_hypergraph_smoke_test.cc`: `AConcurrectHypergraph`
	- \+ `verifyBlockVolumes(hg, k)`, `verifyBlockCutWeights(hg, k)` analog. to `verifyBlockWeightsAndSizes(hg, k)`
	- \+ `verifyConductancePriorityQueue(hg)` (not really..) analog. to `verifyConnectivitySet(hg, k)`
	- `moveAllNodesOfHypergraphRandom(hg, k, obj, show_timing)`: after moving all nodes stats are recomputed (`hypergraph.recomputePartWeights();`) &rArr; call `.recomputePartCutWeights()`, `.recomputePartVolumes()` and `.recomputeConductancePriorityQueue()`
	+ \+ `VerifyBlockVolumesSmokeTest`, `VerifyBlockCutWeightsSmokeTest` analog. to `VerifyBlockWeightsSmokeTest`
	+ \+ `VerifyConductancePriorityQueueSmokeTest` analog. to `VerifyConnectivitySetSmokeTest`

STOPPED HERE: TODO run tests

## Part 2: implementation of the objective function

TODO: write a TODO list for this section :)

### Part 2.0 Conduction priority queue
~~[I'm not sure yet, where to put this pq]~~ &rarr; optional attribute of `PartitionedHypergraph`

**TODO** Mark most methods as *inline* (?)

#### Implementation of PQ
\+ `conductance_pq.h`:
- \+ class `NonnegativeFraction<Numerator, Denominator>` with `operator< , ==, >`, `double_t value()` and getters \+ setters
	- ~~\+ default `operator=(&)`, `operator=(&&)`, `NonnegativeFraction(&&)` to make moving of `ConductancePriorityQueue` &rArr; `PartitionedHypergraph` possible (without `warning`)~~
	- &rarr; default versions are perfect (and are generated as only the trivial constructor is defined)
- `ConductanceFraction := NonnegativeFraction<HyperedgeWeight>`
- \+ class `ConductancePriorityQueue< PartitionedHypergraph > : protected ExclusiveHandleHeap< MaxHeap< PartitionID, ConductanceFraction > >` - addressible max heap with `id = PartitionID`, `key = Conductance`:
	- \+ *private* method `build()`: builds an already filled heap in $\mathcal{O}(k)$
	- \+ `initialize(partitioned_hg, sync)`: initializes underlying `MaxHeap` with `build()` \
		&rarr; \+ `bool initialized`, `bool initialized()` \
		&rarr; \+ `reset(sync)` to return underlying heap to the uninitialized state
	- \+ `globalUpdate(hg, sync)`: a version of `initialize` but for already initialized pq. Should be used after global changes in partitioned hg (e.g. `uncontract(batch, gain_cache)`)
	- \+ *private* metods `lock(sync)`, `uplock(sync)` for a *private* `Spinlock _pq_lock`: \
		Are used by all public methods, when the last parameter `bool synchronized` is set `true`. Per default it is set to `true` only by writing methods except `initialize`, `globalUpdate`, `reset` (These shouldn't be called in parallel).
	- \+ `adjustKey(p, cut_weight, volume, sync)`, `size()`, `bool empty()` - standard pq methods
	- \+ `PartitionID top(sync)`,`PartitionID secondTop(sync)` - return the first and second conductance-wise maximal partitions
	- \+ `vec<PartitionID> topThree(sync)` - returns an **unsorted** vector with 3 top partitions (last elements are `kInvalid`, if `k` < 3). It should help to calculate the gain of a move from $C_i$ to $C_j$ in $\mathcal{O}(3) = \mathcal{O}(1)$ time.
	- \+ `check(phg, sync)` - to check correctness with respect to the phg
	- \+ [potentially useless methods]:
		+ `insert(p, cut_weight, vol, sync)`, `remove(p, sync)`, `deleteTop(sync)` - shouldn't be used if `k = const` 
		+ `topFraction(sync)`, `secondTopFraction(sync)`, `topConductance(sync)`, `secondTopConductance(sync)` - can be emulated witt their `PartitionID` versions + `getCutWeight()`, `getVolume()` from `PartitionedHypergraph`
		+ `getFraction(p, sync)`, `getConductance(p, sync)` shuldn't be used as this information could be obtained from `PartitionedHypergraph`. That way, there will be less syncronization problems withthe underlying pq.

&rarr; `priority_queue.h`: 
- `ExclusiveHandleHeap`:
	- \+ `operator=(&)`, `operator=(&&)`, `ExclusiveHandleHeap(&&)` analog. to the copy constructor of `ExclusiveHandleHeap` \
		&rarr; *Reason:* to make moving of `ConductancePriorityQueue` &rArr; `PartitionedHypergraph` possible (without `warning`)
- `Heap`, `HandleBase`: nothing changed, as only a normal constructor defined &rArr; implicitly defined move assigment operator, move constructor and copy constructor are used.


#### Support of PQ in PartitionedHypergraph

Update of `_conductance_pq` (if enabled): 
1) when partition volume changes
2) when cut weight changes
3) when total volume (potentially) changes


`partitioned_hypergraph.h`:
- \+ `#include "conductance_pq.h"`
- \+ `ConductancePriorityQueue<Self> _conductance_pq` - optional attribute \
	&rarr; \+ `bool hasConductancePQ()`, `enableConductancePQ()`
- \+ `double_t conductance(p)` - calculates conductance of a partition without using `_conductance_pq` \
	returns -1 if volume of partition is 0 &rArr; maybe should throw an exception **???** 
- \+ `bool checkConductancePriorityQueue()` - for testing. True if not initialized.
- \+ `recomputeConductancePriorityQueue()`: runs global update (only for testing (?))

- in 2 constructors: `_conductance_pq()` &lArr; no `initialize()` **!!!**
- `resetData()`: call `_conductance_pq.reset()` in parallel, if enabled
- `uncontract(batch, gain_cache)`: update conductance with `globalUpdate(hg, sync = false)` after `_part_volumes` is finalized (after "// update _part_volumes part 2") \
	**after this the gain cache should be updated (?)** \
    &rArr; `ConductanceGainCache::initializes_gain_cache_entry_after_batch_uncontractions = true` **???**
- `restoreLargeEdge(he)`: call `_conductance_pq.globalUpdate(..)` (running time is already $\Omega(k)$ &rArr; `globalUpdate(..)` is a better choice than `adjustKey(..)` k times)
- `restoreSinglePinAndParallelNets(hes_to_restore, gain_cache)`: call `globalUpdate`: \
	if `_total_volume` is changed (i.e. if a single-pin net is restored), we need to update the whole heap &rArr; $\mathcal{O}(k)$ \
	Else, no update is needed \
	**!!!**: `gain_cache` is updated before `globalUpdate(..)` &rArr; `gain_cache` should have a reference of `_conductance_pq` **???**
- `setNodePart(u, p)` - sets partition for the first time \
	&rArr; `_conductance_pq` shouldn't be initialized yet \
	&rArr; not touched
- `changeNodePart(u, from, to, ...)`: call `adjustKey()` for `from` and `to` \
	!!! update conductance pq after `updatePinCountOfHyperedge(...)` as it updates part cut weight
- `initializePartition()` - not touched for now \
	&rarr; **TODO** initialize `_conductance_pq` somewhere for the case of conductance objective fuction	
- `resetPartition()`: `_conductance_pq.reset()`
- `memoryConsumption`: no info about memory consumption from `priority_queue.h` \
	&rarr; **!!!** for now no info about `ConductancePQ` (TODO **???**)

### Part 2.1 Guide: Setup

1. `partition/context_enum_classes.h`: 
	- \+ `Objective::conductance_local`, `Objective::conductance_global` - new enum types in `Objective`:
		- `conductance_local`: gain of a move from $C_i$ to $C_j$ is the decrease in maximal conductances of the cuts of $C_i$ and $C_j$
		- `conductance_global`: gain of a move from $C_i$ to $C_j$ is the decrease in overall maximal conductances of all cuts (*mostly 0...?*)
	- \+ `GainPolicy::conductance_local`, `GainPolicy::conductance_global` - new enum types in `GainPolicy`:
		analog to `Objective`
2. `partition/context_enum_classes.cpp`: 
	- `operator<< (os, Objective)` and `operator<< (os, GainPolicy)` for `conductance_local`, `conductance_global` (mapping from enum to string)
3. `partition/metrics.cpp`: 

**STOPPED HERE**