# Guide to my code changes
### Notes:
- `HyperedgeWeight` has to be positive (or at least non-negative). Hence, `NonnegativeFraction`. Denominator could be 0: in that case, the `value` is `std::numeric_limits<double_t>::max()`
- `partitioned_hypergraph.h`: `extract(..[4])`, `extractAllBlocks(..[4])`: an extracted hypergraph has no original weighted degrees, original total value **!!!** (could be solved by adding public methods `setNodeOriginalWeightedDegree(u, d)`, `setOriginalTotalVolume(w)` in hypergraphs -- not needed for now, as blocks are used by `recursive_bipartitioning.cpp`)
- `using HypergraphVolume = uint64_t` (defined in `../datastructures/hypergraph_common.h`)

### Problems:
- test `HasCorrectGainsAfterNLevelUncontractionWithLocalizedRefinement` didn't finish in several minutes, but no _conductance_pq.lock(1) was called. Non-stopping runs of actions on github started after implementing original hypergraph stats (original volumes, original weighted degrees, ...)

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
- disable single pin removal for objective conductance (`Objective::conductance_local`,`Objective::conductance_global`) &rArr; check assertions / tests that assure that all single pin nets are removed at coarsening stares etc.

### Next steps:
- implement `_total_volume`, weighted degrees + "their" functions in `static_hypergraph.h/.cpp`: *done*
- write tests (better somewhere in-between): *done*
- implement pq for all conductances (pairs `<Vi_vol, Vi_cut_weight>`)
	- **!!!** Problem: conductance depends on `total_vol` \
	&rarr; heap will not be correct after changes in `total_vol` :( \
	    **???** (= contract, uncontract, remove single-pin?) \
	&rarr; Naive: rebuild heap after changes in `total_vol`? \
	**???** &rarr; implement `increase_total_vol` and refactor to do smth on changes on `total_vol`..?  - *done*
- *STOPPED HERE*
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

### Part 1.1: Hypergraph stats (total volume, weighted degrees)
**???** Also implement for other (heper-)graphs (not only dynamic hg?) as in `lib_generic_impls.h`?
&rarr; for static hg

#### Total Volume:

**For ``DynamicHypergraph``**:

`dynamic_hypergraph.h`:
+ \+ `total_volume` - `std::atomic<HypergraphVolume>` due to changes (analog. to `_contraction_index`)
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
- \+ `_total_volume`: `HypergraphVolume` &rarr; think: better unsigned long long?
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

##### For `DynamicHypergraph`:

&rarr; `dynamic_hypergraph.h`:
- \+ `nodeWeightedDegree(u)` analog. to `nodeDegree()` 
- \+ **!!!** `decreaseNodeWeightedDegree(u, w)`: for dealing with single-pin he (`remove`, `restore`)
- move constructor, move assigment: call `adjustHypergraphPtr(this)` for `_incident_nets` [debug]

&rarr; `incident_net_array.h'
- \+ `CAtomic<HypergraphVolume> weighted_degree` of a node in `Header` of its incident_net's list (analog. to `IncidentNetArray::Header::degree`)
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
    &rarr; for weighted degree, I use `tbb::enumerable_thread_specific< parallel::scalable_vector<HypergraphVolume> >` \
    [analog. to `IncidentNetArray::Header::degree`] \
	**!!!** if `hyperedge_weight_ptr` is passed on from `construct(..)` in `DynamicHypergraphFactory`, we should use it, else use weight=1 for all he instead of `_hypergraph_ptr`, as **the hypergraph is constructed in parallel to its incident net array.** \
	**!!!** Initialize weighted degree to 0 in header(p) as **no `Header` constructor is called** due to pointer tricks with `static_cast` *[debug]*

##### For ``StaticHypergraph``

`static_hypergraph.h`:
- \+ `_weighted_degrees = Array<HypergraphVolume>`
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
		**!!!** use thread-local storage `tbb::enumerable_thread_specific< parallel::scalable_vector < HypergraphVolume> > local_weighted_degree_per_vertex(num_hypernodes, 0);` analog. to `num_incident_nets_per_vertex` [debug]

### Part 1.2: Partitioned hypergraph stats (access to hgInfo, number of pins, volumes, cut weights)

**!!!** Add new methods (except for some testing-only / intern used only) to `partitioned_graph.h` with `UnsupportedOperationException("..")` to make instantiation of `ObjectiveFunction::operator()` possible (templated...)

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
- `memoryConsumption()`: add child for `_part_weights`
- **???** `freeInternalData()` : used by destructor for external memory usages (connectivity...) \
	&rArr; nothing changed


#### Cut weight of Vi:
&rarr; `partitioned_hypergpaph.h`
- \+ `_part_cut_weights` (`vec<CAtomic>` **!!!**)
- constructors, `resetData()` - trivially adjusted analog. to `_part_weights`
- \+ `partCutWeight(p)` - getter
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
- `memoryConsumption(parent)`: add child for `_part_cut_weights` 
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
	- `moveAllNodesOfHypergraphRandom(hg, k, obj, show_timing)`: after moving all nodes stats are recomputed (`hypergraph.recomputePartWeights();`) &rArr; ~~call `.recomputePartCutWeights()`, `.recomputePartVolumes()` and `.recomputeConductancePriorityQueue()`~~ do nothing (or else the new tests `VerifyBlockVolumesSmokeTest`... are useless)
	+ \+ `VerifyBlockVolumesSmokeTest`, `VerifyBlockCutWeightsSmokeTest` analog. to `VerifyBlockWeightsSmokeTest`
	+ \+ `VerifyConductancePriorityQueueSmokeTest` analog. to `VerifyConnectivitySetSmokeTest`

STOPPED HERE: TODO run tests

## Part 2: implementation of the objective function

TODO: write a TODO list for this section :)

### Part 2.0 Conduction priority queue
~~[I'm not sure yet, where to put this pq]~~ &rarr; optional attribute of `PartitionedHypergraph`

**TODO** Mark most methods as *inline* (?)

#### Implementation of The ConductancePriority Queue
##### Nonnegative Fraction
\+ `nonnegative_fraction.h`:
- \+ class `NonnegativeFraction<Numerator, Denominator>` in ``with `operator< , ==, >`, `double_t value()` and getters \+ setters
- ~~\+ default `operator=(&)`, `operator=(&&)`, `NonnegativeFraction(&&)` to make moving of `ConductancePriorityQueue` &rArr; `PartitionedHypergraph` possible (without `warning`)~~
- &rarr; default versions of constructors are perfect (and are generated as only the trivial constructor is defined)
- default value is $0 / 0 = - \infty$
- `denominator` could be 0: in that case, `value() = std::numeric_limits<double_t>::max()`
- comparing cautious, to avoid overflow, division by zero:
	- \+ `bool null()`, `isLess_nullCases(other)`, `isEqual_nullCases(other)` - to cover the case $0 / 0$ \
		**!!!** we say $0 / 0 = -\infty, 0 / x = 0, x / 0 = +\infty$ \ 
		this way, parts with zero weight edges are not disturbed &rarr; TODO: check if suitable 
	- \+ `bool small()`, `isLessQuick(other)`, `isEqualQuick(other)` - by $a/b < c/d \equiv a \cdot d < b \cdot c$ \
		**!!!** only suitable for *small* fractions (i.e. in `uint32_t`, to cause no overflow), not null cases (or else everything is equal to $0/0$ &rArr; heap problems) [debug]
	- \+ `isLess_zeroCases(other)`, `isEqual_zeroCases(other)` - to chech for comparing fractions with zeroes in numerator **or** denumerator
	- \+ `gcd(a, b)`, `reduce(), ``isLessSlow(one, other)`, `isEqualSlow(one, other)` - tricks with integral division and reduction: regulary check for `.._zeroCases(..)`, `small()`

##### Delta Value 
Needed for `changeNodePart(..)` of `partitioned_hypergraph.h`, `adjustKeyByDeltas(..)` of `conductance_pq.h`

\+ `delta_val.h`:
+ \+ `Deltavalue<NonnegativeNumT>` - Delta Value for adjusting keys in the conductance pririty queue during parallel calls of `changeNodePart(u, from, to, ..)`. \ Value `_val` is always non-negative. Boolean flag `_negative` tells the sign.:
	+ \+ getters `bool isNegative()`, `NonnegativeNumT abs()`, `operator<<`
	+ \+ operators `+=, -=, -, +, >, <, ==, !=` for arithmetics with other `DeltaValue` or `NonnegariveNumT`. **!!!** `NonnegativeNumT` shouldn't go first (?)
	- all arithmetical operators use `operator+= (NonnegativeNumT)` or `operator-= (NonnegativeNumT)` which contain an assertion about overflow.

&rarr; `partitioned_hypergraph`: include `delta_val.h` before `conductace_pq.h`

##### ConductancePQ
`hypergraph_commons.h`: [debug: here to avoid include loop]
- include `nonnegative_fraction.h`
- \+ `using ds::ConductanceFraction = NonnegativeFraction<HypergraphVolume>`
- \+ `struct ds::ConductanceInfo { ConductanceFraction, PartitionID}` 

\+ `conductance_pq.h`:
- include `delta_val.h`
- `~~ConductanceFraction := NonnegativeFraction<HypergraphVolume>`~~ moved to `hypergraph_common.h`
- \+ class `ConductancePriorityQueue< PartitionedHypergraph > : protected ExclusiveHandleHeap< MaxHeap< PartitionID, ConductanceFraction > >` - addressible max heap with `id = PartitionID`, `key = Conductance`:
	- \+ *private* method `build()`: builds an already filled heap in $\mathcal{O}(k)$
	- \+ `initialize(partitioned_hg, sync=false)`: initializes underlying `MaxHeap` with `build()` \
		&rarr; \+ `bool initialized`, `bool initialized()` \
		&rarr; \+ `reset(sync=true)` to return underlying heap to the uninitialized state
	- \+ `globalUpdate(hg, sync)`: a version of `initialize` but for already initialized pq. Should be used after global changes in partitioned hg (e.g. `uncontract(batch, gain_cache)`)
	- \+ **public** metods `lock(sync)`, `uplock(sync)` for a *private* `Spinlock _pq_lock`: \
		Are used by all "sync"-versions of public methods, when the last parameter `bool synchronized` is set `true`. Per default it is set to `true` only by writing methods except `initialize`, `globalUpdate`, (These shouldn't be called in parallel). \
		**Not used by normal - const - versions of getters** \
		used by `PartitionedHypergraph::changeNodePart(..)`
	- \+ `size()`, `bool empty()` - standard pq methods
	- \+ `adjustKey(p, cut_weight, volume, sync)` - uses `SuperPQ::adjustKey(..)`.  \ 
		**!!!** `SuperPQ::adjustKey(..)` does nothing when fractions are equal as fractions \
		&rArr; the new numerator and denumerator are not set automaticly, if proportion is the same \
		&rArr; manualy change key in the `SuperPQ::heap` **after** calling `SuperPQ::adjustKey(..)` [debug] \
		**!!!** isn't used by any function (previously used by `changeNodePart`), but due to synchronization problems is replaced by `adjustKeyByDeltas`
	+ \+ ``adjustKeyByDeltas()` - used by `changeNodePart(..)` of `partitioned_hypergraph.h` to accumulate changes in `_part_cut_weights[p]` and `_part_volumes[p]` / `_part_original_volumes[p]` and apply them only once they look 'reasonable' (nonnegative cut weights or part volumes; cut weights or part volumes that are greater than `_total_volume`; cut weights and part volumes that are as sum greater than `_total_volume`) \
	After calling `SuoerPQ::adjustKey(p, f)` we still (as in `adjustLKey(..)`) need to set `SuperPQ::heap[position[p]] = f`
	+ \+ `vec<DeltaV> _delta_part_volumes, _delta_cut_weights` with `using DeltaV = DeltaValue<HypergraphVolume>` - set empty in constructor, filled with `0` in `initialize()`, cleared in `reset()`, asserted as filles with `0` in `flobalUpdate(hg)`, `check(hg)`, `uodateTotalVolume(..)` \
	used for `adjustKeyByDeltas(..)`  to accumulate changes in `_part_cut_weights[p]` and the correct version of part volumes (current or original), that are clearly not yet finished (due to simultaneous calles of `changeNodePart` of `partitioned_hypergraph.h`). This changed are applied once they look *reasonable* (see `adjustKeyByDeltas(..)`)
	- \+ `ConductanceInfo top(sync)`,`ConductanceInfo secondTop(sync)` - return `ConductanceInfo {ConductanceFraction, PartitionID}` of the first, second conductance-wise maximal partitions
	- \+ `vec<ConductanceInfo> topThree(sync)` - returns an **unsorted** vector with `ConductanceInfo {ConductanceFraction, PartitionID}` of the top 3 partitions (last elements are `{0/0, kInvalid}`, if `k` < 3). It should help to calculate the gain of a move from $C_i$ to $C_j$ in $\mathcal{O}(3) = \mathcal{O}(1)$ time.
	- \+ `bool isHeap() const` - version with hushed log
	- \+ `check(phg, sync)` - to check correctness with respect to the phg
	- \+ [potentially useless methods]:
		+ `insert(p, cut_weight, vol, sync)`, `remove(p, sync)`, `deleteTop(sync)` - shouldn't be used if `k = const` 
		+ `topFraction(sync)`, `secondTopFraction(sync)`, `topConductance(sync)`, `secondTopConductance(sync)` - can be emulated witt their `PartitionID` versions + `getCutWeight()`, `getVolume()` from `PartitionedHypergraph`
		+ `getFraction(p, sync)`, `getConductance(p, sync)` shuldn't be used as this information could be obtained from `PartitionedHypergraph`. That way, there will be less syncronization problems withthe underlying pq.
	- \+ `memoryConsumption()`: uses `SuperPQ::memoryConsumption()` and size of `_complement_val_bits`

&rArr; `priority_queue.h`:
+ \+ `ExclusiveHandleHeap<HeapT>::memoryConsumption` - uses `HeapT::size_in_bytes()` of `Heap` (which should be always used as `HeapT`)


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
+ \+ `#include "delta_val.h", "conductance_pq.h"`
+ \+ `ConductancePriorityQueue<Self> _conductance_pq` - attribute

+ \+ `private bool _has_conductance_pq = fasle` - is set to `true` **at the end** of initialization, `_conductance_pq` 
+ \+ `private bool _uses_conductance_pq = true` - if set to `false`, no conductance priority queue is not maintained
+ \+ `bool needsConductancePriorityQueue()` - initializes pq, if it should be done, but wasnt done yet:
	- ~~ **!!!** Called at the end of `partitionImpl()` for all defined partitioners!~~ \
	[undone, because `multilevel.cpp` uses `setOnlyNodePart()` + `initializePartition()` &rArr; `_conductance_pq` is initialized after setting all the `PartitionID`-s]
+ \+ `bool hasConductancePriorityQueue() const` - just returns `_has_conductance_pq`. Is called by all getters to conductance pq (top etc)

+ \+ `bool _conductance_pq_uses_original_stats = true`, `conductancePriorityQueueUsesOriginalStats()` - `enableConductancePriorityQueue()` looks at `_conductance_pq_uses_original_stats` and configures `_conductance_pq`
+ \+ `disableUsageOfOriginalStatsByConductancePriorityQueue()`, `enableUsageOfOriginalStatsByConductancePriorityQueue()`

+ \+ `enableConductancePriorityQueue()` - initializes pq
- \+ `double_t conductance(p)` - calculates conductance of a partition without using `_conductance_pq` \
	returns -1 if volume of partition is 0 &rArr; maybe should throw an exception **???** 
- \+ `bool checkConductancePriorityQueue()` - for testing. True if not initialized.
- \+ `recomputeConductancePriorityQueue()`: runs global update (only for testing (?))
- \+ `ConductanceInfo topPartConductanceInfo()`, `ConductanceInfo secondTopPartConductanceInfo()`, `topThreePartConductanceInfos()` in `PartitionedHypergraph` - initialize `_conductance_pq` if needed (via `needsConductancePriorityQueue()`)
- \+ `conductancePriorityQueue()` to get a const pointer to `_conductance_pq`
- \+ `resetConductancePriorityQueue()`: `_conductance_pq.reset()` + `_has_conductance_pq = false` [debug]

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
- `changeNodePart(u, from, to, ...)`: 
	- call `needcConducatncePriorityQueue()` before the start of changing node part (i.e. if weight was approbed) - so they won't pause updates to initialize pq in the process &rArr; a bit better parallelizm... 
	- call ~~`adjustKey()`~~ `adjustKeyByDeltas(..)` for `from` and `to`  to avoid needing heavy locks and to not break `_conductance_pq` \
		!!! update conductance pq after `updatePinCountOfHyperedge(...)` as it updates part cut weight \
		**!!!** ~~lock `_conductance_pq` when changing part volumes, as `changeNodePart` can be called concurrently &rArr; some threads will be calling `_conductance_pq.adjustKey(..)`, when others are changing part (original & current) volumes [debug] ~~ [debug + adjustKeyByDeltas] 
- &rArr; force `updatePinCountOfHyperedge(..)` to return deltas of `part_cut_weights` pf `from` and `to` as `vec<DeltaValue<hypergraphVolume>>` (`<from, to>`)
- `initializePartition()` - for now initialized pq here--- \
	&rarr; **TODO** initialize `_conductance_pq` somewhere for the case of conductance objective fuction	
- `resetPartition()`: calls `resetConductancePriorityQueue()`
- `memoryConsumption(parent)`: no info about memory consumption from `priority_queue.h` \
	&rarr; **!!!** for now no info about `ConductancePQ` (TODO **???**) \
	&rarr; done with `PriorityQueue::memoryConsumption()`

&rarr; mirrored interface in `partitioned_graph.h`

### External Support of ConductancePQ

in `mt-kahypar/partition/multilevel.cpp`:
- `multilevel_partitioning(..)` : in `## INITIAL PARTITIONING ##` after partitioning call `phg.usesConductancePriorityQueue()` to ensure enabling and initializing of `_conductance_pq` if needed[debug...]

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
	- \+ `ObjectiveFunction<PartitionedHypergraph, Objective::conductance_local>` and `conductance_global` analog.-  template specializations for new `Objective` enum type:
		- &rarr; `operator()(phg, he)` returns 0 if `he` isn't in the most expensive cut, otherwise returns edge weight divided by $min{_part_volume(p), _total_volume - _part_volume(p)}$
		- &rArr; depends on `conductance_pq` of `PartitionedHypergraph`
		- **TODO**: What if several parts have the biggest conductance?
		- **Problem**: value of `ObjectiveFunction` has to be `HyperedgeWeight`: \
		Current solution: `current_multiplier = phg.totalVolume() / phg.k()`. Problems? -> **ASK!!!**
		- if the contribution is too big, a message is printed (`LOG`), the returned value is `std::numeric_limits<HyperedgeWeight>::max()`
	- `contribution(...)`, `quality(...)`: add new objective functions to the switch statements
	+ \+ `HyperedgeWeight compute_conductance_objective(&phg)` - computes conduction without looping through nets
	- in `compute_objective_parallel(..)`, `compute_objective_sequentially(..)` add switch to compute conductance objective without looping through nets 
4. `partition/refinement/gains/gain_definitions.h`:
	- \+ `ConductanceLocalGainTypes`, `ConductanceGlobalGainTypes` - gain type structs analog. to `CutGainTypes`:
		- should in the end contain all relevant implementations for the gain computation in the refinement algorithms &rarr; *to be implemented*
		- for now used `CutGainTypes::[..]`
	- add these classes to the `GainTypes` list
	- and these classes to the macro ~~`INSTANTIATE_CLASS_WITH_TYPE_TRAITS_AND_GAIN_TYPES`~~ [**!!!** no such macro] &rArr; Currently added to all the macros with `CutGainTypes`:
		- `_LIST_HYPERGRAPH_COMBINATIONS`
		- `_INSTANTIATE_CLASS_MACRO_FOR_HYPERGRAPH_COMBINATIONS`
		- `SWITCH_HYPERGRAPH_GAIN_TYPES`
5. `partition/refinement/gains/gain_cache_ptr.h`:
	- `GainCachePtr`: 
		- add new GainPolicy types to all switch statements:
		`applyWithConcreteGainCache(..)`, `applyWithConcreteGainCacheForHG(..)`, `constructGainCache(..)` \
		**!!!** for now used `CutGainCache` &rarr; ***to be changed later**
5. `partition/deep_multilevel.cpp`: 
	- `bipartition_each_block(...)`: ~~add the `GainPolicy` type of new objective functions to the switch statement~~ [no switch statements] &rarr; **nothing changed**
6. `partition/context.cpp`:
	- ~~`sanityCheck(..)`~~ [nothing there] `Context::setupGainPolicy()`: create a mapping between the enum type `Objective` and `GainPolicy` (for `conductance_local` and `clonductance_global`)
7. `partition/registries/register_policies.cpp`: 
	- `== Gain Type Policies ==`: Create a mapping between the enum class `GainPolicy` and its gain type struct (for `conductance_local` and `clonductance_global`)
8. `partition/refinement/gains/bipartitioning_policy.h:`:
	- `useCutNetSplitting(..)` and `nonCutEdgeMultiplier(..)`: add the `GainPolicy` type of new objective functions to the switch statements.	**!!!** **to be rethought later**. ~~For now~~ How it was:
		- `useCutNetSplitting = true`: already cut nets could be cutting nets in the block (if it will be bipartitioned further) &rArr; cannot remove cut nets
		- `nonCutEdgeMultiplier = 1`: otherwise it would change edge weights... Note: in `deep_multilevel.cpp` only `bipartition_each_block(..)` calls `adaptWeightsOfNonCutEdges(..)` and a partitioned hypergraph is built later &rArr; we could recalculate weighted degrees and total volume in the constructor of a partitioned hypergraph **???** &rarr; **no...**: `recursive_bipartitioning.cpp` changes edge weights of a given partitioned hypergraph...  \
	**Rethought**: throw an `UnsupportedOperationException`, as recursive bipartitioning isn't working for conductance
9. \+ `partition/refinement/gains/conductance_local`, `partition/refinement/gains/conductance_global` - folders to that we will later add all relevant gain computation techniques.	

### Part 2.2 ~~Guide:~~ Initial Partitioning

**Problem**: Recursive bipartitioning's invariant 
>  the cut of all bipartitions sum up to the objective value of the initial k-way partition
is not implementable for (scaled) conductance.

&rArr; 2 ways:
1) use singleton partitioning &rArr; change `k` to the number of nodes in the kernel;
2) use recursive bipartitioning with other objective &larr; TODO later (if enough time) - see `partitioner.cpp  partition(...)  #ifdef KAHYPAR_ENABLE_STEINER_TREE_METRIC ...` 

#### Singleton Partitioning:
Add New preset `clustering` with a singleton IP [use commits `a869e6e` "context for cluster & singleton IP set up", `f799400` "singleton IP with k = num nodes of coarsened hg", `04fc118` "fix uncoarsening bug (con info - input num hyperedges)" from https://github.com/adilchhabra/mt-kahypar]:

My changes: 
1) `cluster_preset.ini`:
	```ini
	# main -> refinement -> fm
	...
	r-fm-rollback-parallel=false
	```
2) `presets.cpp`:
	```cpp 
	std::vector<option>	load_clustering_preset() {
		...
		// main -> refinement -> fm
		create_option("r-fm-type", "unconstrained_fm"),
		create_option("r-fm-multitry-rounds", "10"),
		create_option("r-fm-unconstrained-rounds", "8"),
		create_option("r-fm-rollback-parallel", "false"),
		...
	}
	```

##### New clustering preset and new singleton IP

- in `CMakeLists.txt`:
	- \+ option `KAHYPAR_ENABLE_CLUSTERING_FEATURES` \
		[by Adil: `add_compile_definitions`...; by me: `target_compile_definitions`]
	- \+ in `# meta target for library build, which must be built with all features enabled` target compite definitions \+ `KAHYPAR_ENABLE_CLUSTERING_FEATURES` [my idea]
- in `CMakePresets.json` [**my idea**, as no such file by Adil]:
	- enable `KAHYPAR_ENABLE_CLUSTERING_FEATURES` in `default` preset
	- disable `KAHYPAR_ENABLE_CLUSTERING_FEATURES` in `minimal` preset

- in `config/`:
	- \+ `cluster_preset.ini`: 
		- uses `multilevel_coarsener`, only `singleton` IP, no IP refinement, `label_propagation`, `fm` [TODO: shut `fm` down, if `gain_cache` not implemented (?)]
		- [my] set `r-fm-rollback-parallel=false`
	- `large_k_preset.ini`: \+ `i-enabled-ip-algos=0` 
	- all other `.ini` use `initial_partitioning: i-mode=rb` [recursive bipartitioning] &rArr; no changes

- in `mt-kahypar/`:
	- in `partition/`:
		- `context_enum_classes.h`: 
			- \+ `PresetType::cluster`
			- \+ `InitialPartitioningAlgorithm::singleton = 9` (&rArr;  `UNDEFINED = 10`) 
		- `context_enum_classes.cpp`
			- in `operator<< (os, PresetType)`: \+ case `PresetType::cluster` (str: `"cluster"`)
			- in `operator<< (os, mt_kahypar_partition_type_t)`: \+ case `MULTILEVEL_HYPERGRAPH_CLUSTERING` (str: `"multilevel_hypergraph_clustering"`) 
			- in `operator<< (os, InitialPartitioningAlgorithm)`: \+ case `singleton` (str `"singleton"`)
			+ in `presetTypeFromString(string type)`: \+ case `"cluster"`
			+ in `initialPartitioningAlgorithmFromString()`: \+ case `"singleton"`
		- `conversion.cpp`
			- `to_hypergraph_c_type(preset, instance)`: \+ case `PresetType::cluster` [2 times] &rarr; `STATIC_HYPERGRAPH`
			- `to_partition_c_type(preset, instance)`: \+ case `PresetType::cluster` [2 times: for hg and graphs - graphs are my ides ] &rarr; `MULTILEVEL_HYPERGRAPH_PARTITIONING`
		- `partitioner_facade.cpp`:
			- in `check_if_feature_is_enabled(mt_kahypar_partition_type_t type)`: \+ `#ifndef KAHYPAR_ENABLE_CLUSTERING_FEATURES`
			- in lots of functions: \+ `#ifdef KAHYPAR_ENABLE_CLUSTERING_FEATURES: case MULTILEVEL_HYPERGRAPH_CLUSTERING` \
			[in `partition(hg, contex, ..)`, `improve(phg, context, ..)`, `printPartitioningResults(phg, context, ..)`, `serializeCSV(phg, context, ..)`, `serializeResultLine(phg, context, ..)`, `writePartitionFile(phg, filename)]
		- in `initial_partitioning/`: 
			- \+ `singleton_initial_partitioner.h`,  `.cpp`
			- `CMakeLists.txt`: \+ `singleton_initial_partitioner`
		- in `registries/`:
			- `register_initial_partitioning_algorithms.cpp`
				- \+ `#include "../initial_partitioning/singleton_initial_partitioner.h"`
				- \+ define `SingletonPartitionerDispatcher`
				- in `register_initial_partitioning_algorithms()`: \+ register `SingletonPartitionerDispatcher`
			- `register_refinement_algorithms.cpp`:
				- in `getGraphAndGainTypesPolicy(part_type, gain_policy)`: \+ case `MULTILEVEL_HYPERGRAPH_CLUSTERING`
		- in `coarsening/`:
			- `multilevel_uncoarsener.cpp`:
				- in `rebalancingImpl()`: never rebalance by `PresetType::cluster`
	- in `io/`:
		- `command_line_options.cpp`: 
			- add mentioning of preset type cluster: `" - cluster"` 
			- by `"i-enabled-ip-algos"` add singleton IP to the example
		- `presets.cpp`:
			- in `load_large_k_preset()`: by`// main -> initial_partitioning` add entry for `singleton` (`"0"`) 
			- \+ `load_clustering_preset()`: 
				- uses `"multilevel_coarsener"`, `"1" // singleton" IP`, no IP refinement, refinement until no IP improvement by label propagation and fm (but no flows)
				- [my] set `create_option("r-fm-rollback-parallel", false")`
			- in `loadPreset(preset)`: add case `PresetType::cluster` to call `load_clustering_preset()`
		- `partitioning_output.cpp`: [adjust output for `PresetType::cluster`]
			- in `printPartWeightsAndSizes(hg, context)`:
				- \+ `PartitionID num_clusters`: all not-0-weight clusters
				- print `"Num. of clusters = "` for `PresetType::cluster`
				- adjust twice `bool is_imbalanced` to be `false` by `PresetType::cluster`
				- print red `"Number of Imbalanced Blocks = "` only if not `PresetType::cluster`
	- in `utils/`:
		- `cast.h`:
			- in `typeToString(mt_kahypar_partition_type_t)`: add case `MULTILEVEL_HYPERGRAPH_CLUSTERING`
		- `delete.h`:
			- in `delete_partitioned_hypergraph(phg)`: add case `MULTILEVEL_HYPERGRAPH_CLUSTERING` with `ENABLE_CLUSTERING(..)` around it [last: **my idea**] 
	- in `macros.h`:
		- \+ define `ENABLE_CLUSTERING(X) X` if `KAHYPAR_ENABLE_CLUSTERING_FEATURES` defined

- in `include/`:
	- `mtkahypartypes.h`:
		- in `enum mt_kahypar_partition_type_t`: \+ `MULTILEVEL_HYPERGRAPH_CLUSTERING`
		- in `enum mt_kahypar_preset_type_t`: \+ `CLUSTER` - "computes multilevel hypergraph clustering"
	- `lib_helper_functions.h`:
		- `is_compatible(phg, preset)`: \+ case `CLUSTER`	
		- `get_instance_type(phg)`: \+ case `MULTILEVEL_HYPERGRAPH_CLUSTERING`
		- `getget_preset_c_type(preset)`: \+ case `PresetTupe::cluster`
		- `incompatibility_description(phg)`:
			- in case `MULTILEVEL_HYPERGRAPH_PARTITIONING`: \+ compatible with preset `CLUSTER`
			- \+ case `MULTILEVEL_HYPERGRAPH_CLUSTERING`
		- `create_hypergraph(context,...)`: \+ case `PresetType::cluster` [**my idea**, as `mt_kahypar_create_hypergraph(preset, ...)` from `lib/mtkahypar.cpp` was refactored and uses `context` instead of `preset` now] \
		analog.: `create_graph(context, ..)`, `create_partitioned_hypergraph(hg, context, ...)`, : add case `cluster` [**my idea**]
	- `lib_generic_impls.h`:
		- `switch_phg(phg, f)` add case `MULTILEVEL_HYPERGRAPH_CLUSTERING` \ 
		[**my idea**: to inlude `clustering` in refactored `mtkahypar.cpp` `mt_kahypar_write_partition_to_file(partition, ...)`, `mt_kahypar_get_partition(..)`, `mt_kahypar_get_block_weights(..)`, `imbalance(..)`, `mt_kahypar_cut(phg)`, `mt_kahypar_km1(phg)`, `mt_kahypar_soed(phg)`, `mt_kahypar_steiner_tree(..)`]
- in `lib/`:
	- `CMakeLists.txt`: ~~\+ add `target_compile_definitions(mtkahypar PUBLIC KAHYPAR_ENABLE_CLUSTERING_FEATURES)`~~ [not nere now] &rarr; done in `../CMakeLists.txt`
	- `mtkahypar.cpp`:
		- `to_preset_type(preset)`: \+ case `CLUSTER`
		- [my] `from_preset_type(preset)`: \+ case `PresetType::cluster`
		- `mt_kahypar_create_hypergraph(context, ...)`: do nothing here &rarr; add case `cluster` to `create_hypergraph(context, ..)` in `include/lib_helper_functions.h` \
		analog.: `mt_kahypar_create_graph(context, ...)`, `mt_kahypar_create_partitioned_hypergraph(hg, context, ...)`\
		analog. but solved by adding case `MULTILEVEL_HYPERGRAPH_CLUSTERING` to `switch_phg(phg, f)`in `include/lib_generic_impls.h`: \
		`mt_kahypar_write_partition_to_file(partition, ...)`, `mt_kahypar_get_partition(..)`, `mt_kahypar_get_block_weights(..)`, `imbalance(..)`, `mt_kahypar_cut(phg)`, `mt_kahypar_km1(phg)`, `mt_kahypar_soed(phg)`, `mt_kahypar_steiner_tree(..)`
		[**my idea**]

- in `python/`:
	- `CMakeLists.txt`: ~~add `target_compile_definitions(.. KAHYPAR_ENABLE_CLUSTERING_FEATURES)`~~ [`target_link_libraries(mtkahypar_python PRIVATE MtKaHyPar-LibraryBuildSources)` is already used instead &rArr; should be ok, as in global `CMakeLists.txt` a library `MtKaHyPar-LibraryBuildSources` with `KAHYPAR_ENABLE_CLUSTERING_FEATURES` is defined]

##### Setting up k = numNodes() for singleton IP

- in `mt-kahypar/datastructures/`:
	- `partitioned_hypergraph.h`:
		- \+ `setK(k, init_num_hyperedges)`: resets `_part_weights`, `_part_volumes`, `_part_cut_weights`, `part_original_volumes` [debug], `_conductance_pq` [debug], `con_info` [needs `init_num_hyperedges` to reset `con_info`]&rArr; to be called before assigning part_id's
	- `partitioned_graph.h`:
		- \+ `setK(k)`: to be called before assigning part_id's
- in `mt-kahypar/partition/`:
	- `context.cpp`:
		- `setupPartWeights(total_hg_weight)`: \
			if `partition.preset_type == PresetType::cluster` and not `partition.use_individual_part_weights`, set all `context.partition.perfect_balance_part_weights` and `context.partition.max_part_weights` to `std::ceil(total_hypergraph_weight)`
	- `multilevel.cpp`:
		- `multilevel_partitioning(hg, context, ..)`:
			- make `Context& context` argument non-const to be able to change `k` in case of `PresetType::cluster`
			- before `## Coarsening ##` get `HyperedgeID input_he_count = hypergraph.initialNumEdges();`
			- in `## Initial partitioning ##` set `k = phg.initialNumNodes()` in `context` and `phg`
		- in methods of `Multilevel<TypeTraits>` that call `multilevel_partitioning()` (or call methods that call it, etc.): 
			- make `Context& context` argument non-const to be able to call `multilevel_partitioning(hg, context)`, etc. &rArr; change declaration in `.h` \
			[`Multilevel<TypeTraits>::partition(hg, context, ..)`, `Multilevel<TypeTraits>::partition(phg, context, ..)`, `Multilevel<TypeTraits>::partitionVCycle(hg, phg, context, ..)`]
	- `multilevel.h`:
		- make `Context& context` argument non-const to correspond `.cpp` \
			[`partition(hg, context, ..)`, `partition(phg, context, ..)`, `partitionVCycle(hg, phg, context, ..)`]
	- `partitioner.cpp`:
		- `setupContext(hg, context, ..)`: if `PresetType::cluster`, set `context.partition.k = 2` and `context.partition.epsilon = std::numeric_limits<double>::max();` \ 
			[Adil: this determines how the part weights and contraction limits are defined]

**Sanity check**: compiles, passes the test suite

### Part 2.2.1 Side trip: Disabling single-pin net removal
Note: enabling collective sync_update is discussed in the section about Attributed Gain (Problem with `SyncronizedEdgeUpdate`). The new atributes, set-up functions and getters / setters are analogous to single-pin net removal

#### Rationale

**Problem**: removal of single-pin net(s) changes weighted degrees, volumes &rArr; conductances

&rarr; **Idea** *[Adil]*: never remove single-pin nets &rArr; conductance is always right

**Solution structure**: add a member `bool disable_single_pin_net_removal` in `context` and in `StaticHypergraph`, `DynamicHypergraph`

#### ToDo
- new members in `context` (set according to the `Objective`) and hypergraphs (set by `Factory` and internal methods for testing)
- adjust `removeSinglePinAndParallelNets(..)` and ``restoreSinglePinAndParallelNets(..)` (at least an asserting is needed)
- adjust assertions for success of removal (they should be somewhere)
- write (and potentially adjust) tests :(
- debug..

#### Hypergraphs
##### Static Hypergraph
- `mt-kahypar/datastructures`:
	- `static_hypergraph.h`:
		- \+ `void disableSinglePinNetsRemoval()`
		- \+ `bool isSinglePinNetsRemovalDisabled() const`
		- \+ `private bool _disable_single_pin_nets_removal = false`
		- &rArr; adjust ~~copy and~~ move constructor and `operator=` [copy is deleted]
		- `uncontract(..)`: not supported &rArr; do noting :)
	- `static_hypergraph.cpp`:
		- `contract(communities, deterministic)`: [removes single-pin nets from contracted hg] adjust to disable single-pin nets removal: \
			`## Stage 2 ##`: treat single-pin as `size > 1` if single-pin nets removal is disabled
		- `copy()`, `copy(parallel_tag_t)`: copy `_disable_single_pin_nets_removal`
	- `static_hypergraph_factory.h / .cpp`: 
		- `contract(?)`: copy flag
		- `copmactify(SHg)` is not supported &rArr; no changes
##### Dynamic Hypergraph
- `mt-kahypar/datastructures/`:
	- `dynamic_hypergraph.h`:
		- \+ `bool isSinglePinNetsRemovalDisabled() const`
		- \+ `void disableSinglePinNetsRemoval()`
		- \+ `private bool _disable_single_pin_nets_removal = false`
		- &rArr; adjust ~~copy and~~ move constructor and `operator=` [copy is deleted]
	- `dynamic_hypergraph.cpp`:
		- `removeSinglePinAndParallelHyperedges()`: 
			- single-pin nets are not removed if removal disabled
			- but could be removed if are parallel 
		- `restoreSinglePinAndParallelets()`:
			- add an assertion to ensure that only parallel single-pin nets are removed 
			if single-pin nets removal is disabled
		- `copy()`, `copy(parallel_tag_t)`: copy `_disable_single_pin_nets_removal`
	- `hynamic_hypergraph_factory.h / .cpp`: 
		- `compactify(?)`: set flag

##### Partitioned Hypergraph
- `mt-kahypar/datastructures/`:
	- `partitioned_hypergraph.h`:
		- \+ `bool isSinglePinNetsRemovalDisabled() const`
		- `restoreSinglePinAndParallelNets(hes_to_restore, gain_cache)`:
			if single-pin nets removal is disabled, handle restored single-pin nets as parallel nets
		- `extract(block, ..)`: 
			- disable single-pin nets removal for `extracted_block.hg` if is desabled for `_hg`
			- do not remove single-pinn nets (or s-p nets to be) from block!!!
		- `extractAllBlocks(k, ..)`: 
			- analog. to `extract(block, ..)` but in parallel for all blocks
			- do not remove single-pinn nets (or s-p nets to be) from block!!!

#### Graphs
Mirroring interfaces in `static_graph.h`, `dynamic_graph.h`:
- \+ `void disableSinglePinNetsRemoval()`: unsupported
- \+ `bool isSinglePinNetsRemovalDisabled() const`: `false`

Mirroring interfaces in `partitioned_graph.h`:`
- \+ `bool isSinglePinNetsRemovalDisabled() const`: `false`

#### Context
**!!!** At the hypergraph input `hypergraph = io::readInputFile(..)`, single-pin nets are removed if the corresponding argument `bool remove_single_pin_nets = true` is not reset to `false`

&rArr; call `context.setupSinglePinNetsRemoval()` as early as possible after setting `Objective`
- `mt-kahypar/partition`:
	- `context.h`:
		- `CoarseningParameters`: 
			- \+ attribute `bool disable_single_pin_nets_removal = false` 
			- `operator<< (os, CoarseningParameters)`: print out if single-pin nets removal is enabled / disabled
		- `Context`: 
			- \+ `bool disableSinglePinNetsRemoval() const` - getter
			- \+ `void setupSinglePinNetsRemoval()` - to be called after Objective is set up
	- `partitioner.cpp`:
		- `setupContext(&hg, &context, ..)`: 
			- before all other: if conductance `Objective`, disable single-pin nets in `context`
			- after setup calls [as `context.partition.instance_type` could be initializes there]: if `hypergaph` (*not graph*), disable single-pin nets in `hg`
- `mt-kahypar/application/`:
	- `mt_kahypar.cc`:
		- `main(argc, argv)`:
			- in `hypergraph = io::readInputFile(..)`: set option `bool remove_single_pin_nets = !context.coarsening.disable_single_pin_nets_removal` (is `true` by default)
			- before that call `context.setupSinglePinNetsRemoval()` 
- `mt-kahypar/io/`:
	- `command_line_options.cpp`: 
		- in `processCommandLineInput(&context, ..)` call `context.setupsetupSinglePinNetsRemoval()` after `context.partition.objective = objectiveFromString(s);`

#### Assertions about edge size
- `mt-kahypar/partition/coarsening/multilevel_vertex_pair_rater.h`: - used in by `multilevel_coarsener.h`
	- `fillRatingMap(hg, u, tmp_ratings)`: `ASSERT(edge_size > 1)`:
		- `ScorePolicy::score(edge_weight, edge_size)` needs `edge_size != 1` &rArr; adjust assertion and skip single-pin nets

#### Serialization of hg, context
`/mt-kahypar/io/sql_plottools_serializer.cpp`:
- \+ write out `context.coarsening.disable_single_pin_nets_removal` to serialize `context` correctly \
	[debug: needed in `mt-kahypar/tests/io/sql_plottools_serializer_test.cc` `ASqlPlotSerializerTest.ChecksIfSomeParametersFromContextAreMissing`, removed whitespaces in empty line in `context.h` - they break serialization!]

Also done:
- in `mt-kahypar/partition/context.cpp`:
	- `sanityCheck(..)`: if conductance `Objective`, but not `PresetType::cluster`, throw an `UncupportedOperationException`. 

- *Call hierarchy*:
		- `Context::sanityCheck(..)`
		- &rarr; `partitioner.cpp setupContext(hg, context, ..)`
		- &rarr; `partitioner.cpp Partitioner<TypeTraits>::partition(hg, context, ..)`
		- &rarr; `partitioner_facade.cpp internal::partition(hg, context, ..)`
		- &rarr; `partitioner_facade.cpp partition(hg, context, ..)`
		- &rarr; `lib_helper_functions.h lib::partition_impl(hg, context, ..)`
		- ... &rArr; should be the right place to set `disable_single_pin_nets_removal`..


#### Notes:
- `mt-kahypar/partition/refinement/fm/global_rollback.cpp`: 
	- `recalculateGainForGraphEdgeViaAttributedGains(phg, FMSharedData, e)` ignores single-pin net `e`

#### Tests 
**ToDo**:
- Dynamic, static hg: remove single <..> with disabled removal, copy test
- static: contract
- partitioned: extract block, all blocks
- + run a smoke test with disabled SPNR to catch failed assertions

- `static_hypergraph_test.cc`:
	- \+ `ComparesSinglePinNetsRemovalOptionIfCopiedParallel1 , 2` and `ComparesSinglePinNetsRemovalOptionIfCopiedSequential1 , 2`: tests `desable` + `copy`
	- \+ `ContractsCommunitiesWithDisabledSinglePinNetRemoval1 .. 3`: tests `desable` + `contract`

- `dynamic_hypergraph_test.cc`:
	- \+ `ComparesSinglePinNetsRemovalOptionIfCopiedParallel1 , 2` and `ComparesSinglePinNetsRemovalOptionIfCopiedSequential1 , 2`: tests `desable` + `copy`
	- \+ `RemovesOnlyParallelNets1 , 2`: tests `disable` + `removeSPNandPN`
	- \+ `RestoresOnlyParallelNets1 , 2`: tests `disable` + `restoreSPNandPN`

- `dynamic_partitioned_hypergraph_test.cc`:
	- \+ `ComputesPinCountsCorrectlyIfWeRestoreOnlyParallelNets`: `disable` + `removeSPNandPN` + `initializeGainCache` + `restorePNandPN`

- `partitioned_hypergraph_test.cc`:
	- \+ `ExtractBlockZeroWithCutNetSplittingAndSinglePinNets`, `ExtractBlockOneWithCutNetSplittingAndSinglePinNets`, `ExtractBlockTwoWithCutNetSplittingAndSinglePinNets`: `extract` + `disable`
	- \+ `ExtractAllBlockBlocksWithCutNetSplittingAndSinlePinNets`: `extractAll` + `disable`

### Part 2.2.2 Side trip: Problem: volumes are still not preserved
Contraction of 2 nodes decreases volume by the sum of weights of their shared nets

#### Solution:
1) introduce `_original_weighted_degrees` of a node and `_original_total_volume` to `static_hypergraph.h / .cpp`, `dynamic_hypergraph.h / .cpp`:
	- not changed by removal / restoration of a single-pin net;
	- summed by contraction, subtracted by uncontraction;
2) introduce `_part_original_volumes` to `partitioned_hypergraph.h`:
	- not changed by restoration of single-pin net;
	- ~~substracted (?) by uncontractiion~~
	- adjust by `changeNodePartition(..)`, `setNodePartition(..)`
3) adjust `conductance_pq.h / .cpp`:
	- use per default `_part_original_volumes`, `_original_total_volume`
	- add assertion by `updateTotalVolume(..)`, `globalUpdate` (?)
	- leave a possibility to use current `_part_volumes` and `_total_volume`
4) tests...

#### Original Weighted Degree and Original Total Volume in Hypergraph
##### StaticHypergraph
- `static_hypergraph.h`:
	+ \+ `Array<HypergraphVolume> _original_weighted_degrees`
	+ \+ `nodeOriginalWeightedDegree(u)` - getter
	+ \+ `HypergraphVolume _original_total_volume`
	+ \+ `otiginalTotalVolume()` - getter
	- adjust constructor, move constructor and move assigment operator
	- in `removeEdge(he)`, `removeLargeEdge(he)`, `restoreEdge(he)` and `restoreLargeEdge(he)`: **not** adjust `_original_weighted_degree[pin]` and `_original_total_volume` (used only in tests)
- `static_hypergraph.cpp`:
	- in `contarct(communities, ..)`, stage 4: 
		1) resize `hypergraph._original_weighted_degrees` analog. to `hypergraph._weighted_degrees`
		2) accumulate `hypergraph._original_weighted_degrees` as sum of original weighted degrees of the contracted nodes
		3) `hypergraph._original_total_volume = _original_total_volume`
	- in `copy(parallel_tag_t)`, `copy()`: 
		- copy `_original_total_volume`
		- copy `_original_weighted_degrees` analog. to `_weighted_degrees` [in parallel version: in parallel]
	- in `memoryConsumption(parent)`: analog. to `_weighted_degrees`, add child `"Original Weighted Degrees"` 
- `static_hypergraph_factory.cpp`:
	- in `construct(..)`:
		1) resize `_original_weighted_degrees` at the beginning
		2) after computing weighted degrees, copy it in parallel to `_original_weighted_degrees`
		3) set `_original_total_volume` after calling `computeAndSetTotalVolume()`	

##### DynamicHypergraph
- `incident_net_array.h`:
	- \+ `HypergraphVolume Header::original_weighted_degree` - not atomic!!!
	- \+ `nodeOriginalWeightedDegree()`
	- \+ `setOriginalWeightedDegree()` - for `DynamicHypergaph::compactify(..)`
- `incident_net_array.cpp`:
	- `contract(u, v, ..[3])`: sum up original weighted degrees of contracted nodes
	- `uncontract(u, v, ..[4])`: subtract the original weighted degree of returning hn `v` from the original weighted degree of `u`
	- `construct(edge_vector, hyperedge_weight_ptr)`: assign `original_weighted_degrees` in parallel after calculating `weighted_degrees`
	- **no more changes here!!!**

- `dynamic_hypergraph.h`:
	+ \+ `HypergraphVolume originalTotalVolume` - not atomic as constant
	+ \+ `originalTotalVolume()` - getter
	+ ~~\+ `setOriginalTotalVolume(original_total_volume)` - setter: **!!!** only for `DynamicHypergraphFactory::compactufy(hypergraph)`~~ `DynamicHypergraphFactory` is a *friend* class &rarr; no need :)
	+ \+ `nodeOriginalWeightedDegree(u)`, `setNodeOriginalWeightedDegree(u, w)` - calling the corresponding methods of `InsidentNetArray`
	- in constructor, move constructor, move assigment operator set `_original_total_volume`
- `dynamic_hypergraph.cpp`:
	- `copy(parallel_tag_t)`, `copy()`: copy `_original_total_volume`
- `dynamic_hypergraph_factory.cpp`:
	- `construct(..[5])`: assign `original_total_volume` after calling `updateTotalVolume(parallel_tag_t())`
	- `compactify(hypergraph)`: 
		- `setNodeOriginalWeightedDegree(..)` analog. to `setComunityID(..)` after constructing `compactified_hypergraph`
		- set the `_original_total_volume`

##### Partitioned hypergraph
- `partitioned_hypergraph.h`:
	+ \+ `nodeOriginalWeightedDegree(u)`
	+ \+ `originalTotalVolume()`

#### Original Part Volumes
##### Partitioned Hypergraph
- `partitioned_hypergraph.h`:
	+ \+ `vec< CAtomic<HypergraphVolume> > _part_original_volumes` - volume according to the original weighted degrees for all blocks
	+ \+ `partOriginalVolume(p)`
	+ \+ `recomputePartOriginalVolume()` - uses weighted degrees, method is used for testing
	+ \+ `incrementOriginalVolumeOfBlock(p, w)`, `decrementOriginalVolumeOfBlock(p, w)` - analog. to `incrementVolumeOfBlock(p, w)`...
	+ \+ `initializeBlockOriginalVolumes()`
	+ \+ `applyPartOriginalVolumeUpdates(vec<HypergraphVolume>&)` - for now, needed only for `initializeBlockOriginalVolumes()`
	+ \+ `double_t originalConductance(p)`
	+ \+ `bool _conductance_pq_uses_original_stats = true`, `conductancePriorityQueueUsesOriginalStats()` - `enableConductancePriorityQueue()` looks at `_conductance_pq_uses_original_stats` and configures `_conductance_pq` accordingly
	+ \+ `disableUsageOfOriginalStatsByConductancePriorityQueue()`, `enableUsageOfOriginalStatsByConductancePriorityQueue()`
	- set in 2 constructors, reset in `resetData()` 
	- `setK(k, ..[1])` - resize `_part_original_volumes` !!! [debug]
	- `uncontract(batch, gain_cache)`, `restoreLargeEdge(he)`, `restoreSinglePinAndParallelNets(..)`: 
		- no update of `_part_original_volume` needed
		- if **not** `conductancePriorityQueueUsesOriginalStats()`: run `conductance_pq.globalUpdate()`
	- `setnodePart(u, p)`: increment`_part_original_volumes[p]` by `originalWeightedDegree(u)` by calling the corresponding method
	- `changeNodePart(u, from, to, ..)`: adjust `_part_original_volumes[from \ to]` analog. to `_part_volumes`, adjust keys in `_conductance_pq` according to `conductancePriorityQueueUsesOriginalStats()`
	- `initializePartition()`: call `initializeBlockOriginalVolumes()` in parallel
	- `resetPartition()`: reset `_part_original_volumes`
	- `memoryConsumption(parent)`: add child for `_part_original_volumes` and volumes - TODO
	- `extract(..[4])`, `extractAllBlocks(..[4])`: no changes for now &rArr; an extracted hypergraph has no original weighted degrees, original total value **!!!** (could be solved by adding public methods `setNodeOriginalWeightedDegree(u, d)`, `setOriginalTotalVolume(w)` in hypergraphs -- not needed for now, as blocks are used by `recursive_bipartitioning.cpp`)

##### Conductance Priority Queue
- `conductance_pq.h`:
	+ \+ `_uses_original_stats = true` per default
	+ \+ `disableUsageOfOriginalHGStats()`, `enableUsageOfOriginalHGStats()`
	+ \+ `private getHGTotalVolume(&hg)`, `private getHGPartVolume(&hg, p)` to get original / current stats from `hg`
	+ \+ `private getHGPartCutWeight(&hg, p)` - for uniform communication with the hypergraph
	
	- `initialize(hg, sync)`: use `getHG..` to get total volume, part volumes, part cut weights
	- `reset(sync)`: **doesn't** reset `_uses_original_stats` to `true`
	- `globalUpdate(sync)`: 
		- add an assertion: if uses original stats, original total volume shouldn't have changed
		- use `getHG..` to get total volume, part volumes, part cut weights
	- `check(hg)`: use `getHG..` to get total volume, part volumes, part cut weights
	- `updateTotalVolume(new_total_volume, sync)`: add an assertion, that original stats are not used

#### Tests
- `static_hypergraph_test.cc`:
	- `HasCorrectStats`: check `originalTotalVolume()`
	+ \+ `VerifiesVertexOriginalWeightedDegrees()`: at the begining the same as weighted degrees [analog. to `VerifiesVertexWeightedDegrees`]
	- `RemovesAHyperedgeFromTheHypergraph1 .. 4`: checks unchanged original weighted degrees, original total volume
	- `ComparesStatsIfCopiedParallel`, `ComparesStatsIfCopiedSequential`: checks unchanged original total volume
	- `ComparesWeightedDegreesIfCopiedParallel`, `ComparesWeightedDegreesIfCopiedSequential`: checks unchanged original weighted degrees
	- `ContractsCommunities1 .. 3` : checks unchanged original total volume, recalculated original weighted degrees
	- `ContractsCommunitiesWithDisabledSinglePinNetRemoval1 .. 3`: checks unchanged original total volume, recalculated original weighted degrees
	- `ContractsCommunitiesWithDisabledHypernodes`, `ContractsCommunitiesWithDisabledHyperedges`: checks **unchanged** original total volume (nodes are disabled in contraction, edges - in removal of parallel / single pin nets &rArr; original total volume shouldn't be changed) and recalculated original weighted degrees (with no *accomodation* for disabled nodes / edges...) 

- `dynamic_hypergraph_test.cc`:
	- `HasCorrectStats`: check original total volume
	+ \+ `VerifiesVertexOriginalWeightedDegrees` analog. to `VerifiesVertexWeightedDegrees`
	- `ComparesStatsIfCopiedParallel`, `ComparesStatsIfCopiedSequential`, : check unchanged original total volume
	- `PerformsAContraction1 .. 5`, `PerformAContractionsInParallel1 .. 3`: check unchanged original total volume, summed original weighted degrees
	- `verifyEqualityOfDynamicHypergraphs(..)`: check equality of original total volumes and original weighted degrees
	- `RemovesSinglePinAndParallelNets1 .. 2`, `RemovesOnlyParallelNets1`: check unchanged original total volumes and the same way recalculated original weighted degrees
	- `RestoreSinglePinAndParallelNets1 .. 2`, `RestoreOnlyParallelNets1 .. 2`: check unchanged original total volumes and the same way recalculated original weighted degrees
	- `GeneratesACompactifiedHypergraph1 ..`: check unchanged original total volumes and the same way recalculated original weighted degrees

- `partitioned_hypergraph_test.cc`:
	+ \+ `HasCorrectPartOriginalVolumes`, `HasCorrectPartOriginalVolumesIfOnlyOneThreadPerformsModifications`
	- `PerformsConcurrentMovesWhereAllSucceed`: check unchanged original part original volumes 
	+ \+ `ChecksConductancePQWithOriginalStatsAfterConcurrentMoves`, `ChecksConductancePQWithCurrentStatsAfterConcurrentMoves`
	- `ComputesPartInfoCorrectIfNodePartsAreSetOnly`: check part original volumes 

- `dynamic_partitioned_hypergraph_test.cc`:
	- no changed tests

- `partitioned_hypergraph_smoke_test.cc`:
	+ \+ `verifyBlockOriginalVolumes`: analog. to `verifyBlockVolumes`
	+ \+ `VerifyBlockOriginalVolumesSmokeTest`: analog. to `VerifyBlockVolumesSmokeTest`
	- `VerifyConductancePriorityQueueSmokeTest`: 
		- **rename** to `VerifyConductancePriorityQueueWithOriginalStatsSmokeTest`
		- add assertion about using original stats
	+ \+ `VerifyConductancePriorityQueueWithCurrentStatsSmokeTest`

### Part 2.3 Guide: Label Propagation Refinement

Label propagation algorithm iterates over all nodes in parallel and moves each node to the block **with the highest gain**. 

The algorithm requires two gain techniques: 
1) A *gain computation* algorithm to compute **the highest gain move for a node** 
2) the *attributed gain* technique to **double-check the gain of a node move at the time performed on the partition**. 

#### Guide: outer changes

`partition/`:
- `refinement/`:
	- `gains/gain_definitions.h`:
		+ include `conductance_local_gain_computation.h`, `conductance_global_gain_computation.h`, ~~`conductance_local_attributed_gain.h` and~~ `conductance_global_attributed_gain.h`
		+ replace `GainComputation` and `AttributedGains` members of `ConductanceLocalGainTypes` and `ConductanceGlobalGainTypes`. **!!!** both use `ConductanceGlobalAttributedGains`

#### Synchronized Update (my changes)
##### Intro from guide

- The ```SynchronizedEdgeUpdate``` structs contains the following members:
```cpp
struct SynchronizedEdgeUpdate {
  HyperedgeID he;
  PartitionID from;
  PartitionID to;
  HyperedgeID edge_weight;
  HypernodeID edge_size;
  HypernodeID pin_count_in_from_part_after;
  HypernodeID pin_count_in_to_part_after;
  PartitionID block_of_other_node; // only set in graph partitioning mode
  mutable ds::Bitset* connectivity_set_after; // only set when optimizing the Steiner tree metric
  mutable ds::PinCountSnapshot* pin_counts_after; // only set when optimizing the Steiner tree metric
  const TargetGraph* target_graph; // only set when optimizing the Steiner tree metric
  ds::Array<SpinLock>* edge_locks; // only set when optimizing the Steiner tree metric
};
```
When we move a node from its *source* (```from```) to a *target* block (```to```), we iterate over all hyperedges, perform syncronized data structure updates and call `gain` function for each incident hyperedge of the moved node. The sum of all calls to this function is the attributed gain of the node move.

###### Changes
&rArr; `hypergraph_common.h`:
- `SynchronizedEdgeUpdate`: \+ New attributes  be able to **estimate** AttributedGain (never exact due to rounding. scaling...):
	0) `PartitionID k`
	1) `HypergraphVolume cut_weight_from_after, cut_weight_to_after` 
	2) `vec<ds::ConductanceInfo> top_three_conductance_info_before` - ConductanceInfo of the top 3 conductance-wise parts before the move: [0, 1/2, 2/1]
	3) `HypergraphVolume volume_from_after, volume_to_after`
	5) `HypergraphVolume weighted_degree` (the used version!!!)
	6) `HypergraphVolume total_volume` (the used version!!!)
	7) `PartitionID from, to` 
- \+ f~~orward declaration~~ definition of `ds::ConductanceInfo` [debug: here to avoid include loop]

Initialize new members of `SynchronizedEdgeUpdate` (looked for mentions: skipped lambdas, tests):
- `partitioned_hypergraph.h`:
	- `changeNodePart(u, from, to, ..)`: 
		- after `// Construct sync_update` (before iterating through incident nets) set **if `needsConductancePriorityQueue()`** :
			- `k`, `top_three_conductance_info_before`
			-  `volume_from_after`, `volume_to_after`, `weighted_degree`, `total_volume` [used version only]
		- after iterating through incident nets 
			- set `cut_weight_to_after`, `cut_weight_from_after`
		- if collective sync_update is enabled:
			- call `notify_func(sync_update)` before loop over incident edges
			- call `delta_func(sync_update)` after loop + setting the rest of `sync_update`

- `global_rollback.cpp`:
	- `recalculateGainForHyperedgeViaAttributedGains(&phg, FMShareData&, &e)`, `recalculateGainForGraphEdgeViaAttributedGains()`: 
		- **do nothing**, as used by `recalculateGains(..)` &larr; `revertToBestPrefixParallel(..)`, but I can't use it (parallel computation of attributed gain is incorrect...) 

#### Attributed Gain
##### Intro from guide
The gain of a node move can change between its initial calculation and execution due to concurrent node moves in its neighborhood. 

&rArr; We double-check the gain of a node move at the time performed on the partition via synchronized data structure updates. This technique is called *attributed gains*:
- The *label propagation algorithm* **reverts node moves that worsen the solution quality** by checking the attributed gain value. 
- The attributed gain function implements the following interface:
	```cpp
	static HyperedgeWeight gain(const SynchronizedEdgeUpdate& sync_update);
	```
##### Problem: I can only calculate collective AttributedGain of a move (not a sum for all edges)
&rarr; analog. to `_disable_single_pin_nets_removal`

###### Hypergraphs
####### Static Hypergraph
- `mt-kahypar/datastructures`:
	- `static_hypergraph.h`:
		- \+ `void enableCollectiveSyncUpdates()`
		- \+ `bool areCollectiveSyncUpdatesEnabled() const`
		- \+ `private bool _enable_collective_sync_update = false`
		- &rArr; adjust ~~copy and~~ move constructor and `operator=` [copy is deleted]
	- `static_hypergraph.cpp`:
		- `contract(communities, deterministic)`: set flag of the constructed hypergraph
		- `copy()`, `copy(parallel_tag_t)`: copy `_enable_collective_sync_update`
	- `static_hypergraph_factory.h / .cpp`: `copmactify(SHg)` is not supported &rArr; no changes
####### Dynamic Hypergraph
- `mt-kahypar/datastructures/`:
	- `dynamic_hypergraph.h`:
		- \+ `bool areCollectiveSyncUpdatesEnabled() const`
		- \+ `void enableCollectiveSyncUpdates()`
		- \+ `private bool _enable_collective_sync_update = false`
		- &rArr; adjust ~~copy and~~ move constructor and `operator=` [copy is deleted]
	- `dynamic_hypergraph.cpp`:
		- `copy()`, `copy(parallel_tag_t)`: copy `_enable_collective_sync_update`
	- `hynamic_hypergraph_factory.h / .cpp`: 
		- `compactify(..)`: set flag to the constructed hypergraph

####### Partitioned Hypergraph
- `mt-kahypar/datastructures/`:
	- `partitioned_hypergraph.h`:
		- \+ `bool collectiveSyncUpdatesEnabled() const`, `void enableCollectiveSyncUpdates()`:			
			- `changeNodePart(..)`: if collective sync_update is enabled
				- call `notify_func(sync_update)` in before changing pin counts
				- call `delta_func(sync_update)` after changing pin counts
			- `updatePinCountOfHyperedge(..)`: call `notify_func(sync_update)`, `delta_func(sync_update)` only is collective sync_update is disabled
		- `extract(block, ..)`: 
			- set the flag for `extracted_block.hg` if is desabled for `_hg`
		- `extractAllBlocks(k, ..)`: 
			- analog. to `extract(block, ..)` but in parallel for all blocks
		+ add assertions to `topThreePartConductanceInfos()`

###### Graphs
Mirroring interfaces in `static_graph.h`, `dynamic_graph.h`:
- \+ `void enableCollectiveSyncUpdates()`: unsupported
- \+ `bool areCollectiveSyncUpdatesEnabled() const`: `false`

Mirroring interfaces in `partitioned_graph.h`:`
- \+ `bool collectiveSyncUpdatesEnabled() const`: `false`
- \+ `void enableCollectiveSyncUpdates()`: unsupported

###### Context

call `context.setupCollectiveSyncUpdates()` as early as possible after setting `Objective`
- `mt-kahypar/partition`:
	- `context.h, cpp`:
		- `PartitioningParameters`: 
			- \+ attribute `bool enable_collective_sync_updates = false` 
			- `operator<< (os, PartitioningParameters)`: print out if single-pin nets removal is enabled / disabled
		- `Context`: 
			- \+ `bool enableCollectiveSyncUpdates() const` - getter
			- \+ `void setupSyncUpdatePreference()`: if conductance objective, enable collective sync_update (if not yet) and write in `LOG` [analog. to `setupSinglePinNetsRemoval()`]
			- `sanityCheck()` also sets the flag in context
	- `partitioner.cpp`:
		- `setupContext(&hg, &context, ..)`: 
			- before all other: if conductance `Objective`, set up sync_update preference in `context`
			- after setup calls [as `context.partition.instance_type` could be initializes there]: if `hypergaph` (*not graph*), enable sync_updates in `hg` (if context says so)

- `mt-kahypar/io/`:
	- `command_line_options.cpp`: 
		- in `processCommandLineInput(&context, ..)` call `context.setupSyncUpdatePreference()` after `context.partition.objective = objectiveFromString(s);`

###### Serialization of hg, context
`/mt-kahypar/io/sql_plottools_serializer.cpp`:
- \+ write out `context.partition.enable_collective_sync_updates` to serialize `context` correctly \
	[debug: needed in `mt-kahypar/tests/io/sql_plottools_serializer_test.cc` `ASqlPlotSerializerTest.ChecksIfSomeParametersFromContextAreMissing`, removed whitespaces in empty line in `context.h` - they break serialization!]

##### Conductance Global
**As far as I know** [*to be disproven by failing quality assertions...*], `contribution(phg, he)` is called only by `partitioner.cpp` in the `POSTPROCESSING` phase: `restoreLargeHyperedges()`. So `AttributedGain` is compared to the `quaality(phg)`, which can be (and now is) calculated almost exact (rounded max conductance). Therefore to calculate `AttributedGains` for conductance (both use global version), I compute the difference on new and old rounded conductances (and do this by enabled **collective** sync_updates)

\+ `conductance_global_attributed_gain.h`:
+ \+ `ConductanceGlobalAttributedGains`:
	+ \+ `public HyperedgeWeight compute_conductance_objective(..)` - computes scaled + rounded conduction exactly the same way as `quality(const &phg)` in `metrix.cpp` 
	- `static HyperedgeWeight gain(&sync_update)` - computes new fractions for `to` and `from`, finds the biggest other fraction (in `sync_update.top_three_conductance_info_before` - vector with $\le 3$ biggest conductance fractions), computes old and new top conductances via `HyperedgeWeight compute_conductance_objective(..)`, returns `new_conductance - old_conductance`

##### Conductance Local
Uses the same attributed gains
&rArr; in `gain_definitions.h` set `ConductanceLocalGainTypes::AttributedGains = ConductanceGlobalAttributedGains`

~~\+ `conductance_local_attributed_gain.h`:~~
~~+ \+ `ConductanceLocalAttributedGains`~~


#### Gain Computation
##### Intro From Guide

All gain computation techniques inherit from ```GainComputationBase``` (see ```partition/refinement/gains/gain_compute_base.h```). The base class has two template parameters: the derived class and the attributed gain implementation of the objective function (curiously recurring template pattern, avoids vtable lookups). The base class calls the ```precomputeGains(...)``` function of the derived class, which has the following interface:
```cpp
template<typename PartitionedHypergraph>
void precomputeGains(const PartitionedHypergraph& phg,
                     const HypernodeID hn,
                     RatingMap& tmp_scores,
                     Gain& isolated_block_gain);
```
We split the gain computation in two steps: 
1) compute the gain of moving the node into an isolated *= non-adjasent* block (stored in ```isolated_block_gain```) [0 by me, as they are different for different non-adjasent blocks &rArr; have to be in `tmp_scores` - analog. to steiner tree]
2) moving the node from the isolated block to all adjacent blocks (stored in ```tmp_scores```). ```tmp_scores``` can be used similar to an ```std::vector``` and has exactly k entries. The gain of moving a node to a block ```to``` can then be computed by ```isolated_block_gain - tmp_scores[to]``` [**not for conductance: there just `tmp_scores[to]`**] (a negative value means that moving the node to block ```to``` improves the objective function). However, the derived class still implements a function that computes the gain to a particular block:
```cpp
HyperedgeWeight gain(const Gain to_score, // tmp_scores[to]
                     const Gain isolated_block_gain);
```
The implementation of this function is most likely ```isolated_block_gain - to_score``` (except for the Steiner tree metric, where it is `to_score` &rArr; **by conductance, too**).

##### Conductance Global
All lookups for gain of a node in `gain_computation_base.h` pass through `gain(to_score, isolated_block_gain)` &rArr; should be ok to set `isolated_block_gain = 0`

\+ `conductance_global_gain_computation.h`
+ \+ `ConductanceGlobalGainComputation`:
	+ `private ds::Bitset _all_blocks`, `private tbb::enumerable_thread_specific<ds::Bitset> _local_adjacent_blocks` - bit-sets for iterating over concidered blocks [took from `steiner_tree_gain_computation.h`], set in constructor, changed in `changeNumberOfBlocksImpl(new_k)`
	+ `static constexpr size_t BITS_PER_BLOCK = ds::StaticBitset::BITS_PER_BLOCK` - also took from `steiner_tree_gain_computation.h`
	+ \+ [my] `setPreSyncUpdataData<artitionedHypergraph>(phg, ...)` - sets `to`, `volume_to_after`, `cut_weight_to_after`, `cut_weight_from_after` - almost all the data needed for the emulation of sync_update. Note: cut weights have to be adjusted later in the `precomputeGains(..)` method
	- `void precomputeGains(const &phg, const hn, &tmp_scores, &isolatedBlockGain, consider_non_adjacent_blocks)`: 
		- analog. to steiner tree gets a "list" of all considered blocks;
		- gets current conductance objective value via `compute_conductance_objective(..)` of `ConductanceGlobalAttributedGains`
		- for each considered block:
			- collects data needed for simulation of `sync_update` for a conductance objective (via `void setPreSyncUpdataData(&..)` + cut weights adjusted in a loop);
			- calls `ConductanceGlobalAttributedGains::gain(sync_update)` to get gain
			- sets `tmp_scores[to] = conductance_after - conductance_now`
	- `HyperedgeWeight gain(to_score, isolated_block_gain)`: returns `to_score`

##### Conductance Local
Computes the "gain" of move as change in maximal conductance between to and from blocks. Otherwise analog. to Conductance Global

\+ `conductance_local_gain_computation.h`:
- \+ `ConductanceLocalGainComputation`
	+ `private ds::Bitset _all_blocks`, `private tbb::enumerable_thread_specific<ds::Bitset> _local_adjacent_blocks` - bit-sets for iterating over concidered blocks [took from `steiner_tree_gain_computation.h`], set in constructor, changed in `changeNumberOfBlocksImpl(new_k)`
	+ `static constexpr size_t BITS_PER_BLOCK = ds::StaticBitset::BITS_PER_BLOCK` - also took from `steiner_tree_gain_computation.h`
	+ \+ [my] `setPreIncidentEdgeLoopData<artitionedHypergraph>(phg, ...)` - sets `cut_weight_from_before`, `cut_weight_to_before`, `cut_weight_from_after`, `cut_weight_to_after`, `volume_to_before` [used version], `volume_to_after`[used version]: \
	&larr; almost all the data needed for the emulation of a move. \
	Note: cut weights have to be adjusted later in the `precomputeGains(..)` method
	- `void precomputeGains(const &phg, const hn, &tmp_scores, &isolatedBlockGain, consider_non_adjacent_blocks)`: 
		- analog. to steiner tree gets a "list" of all considered blocks;
		- gets current from fraction
		- for each considered block:
			- gets current to fraction 
			- &rarr; computes current local conductance objective value of maximal fraction between from and to via `compute_conductance_objective(..)` of `ConductanceGlobalAttributedGains`  (at least one part isn't empty &rArr; max. isn't $0 / 0$)
			- collects data needed for simulation for move (via `void setPreIncidentEdgeLoopData(&..)` + cut weights adjusted in a loop);
			- gets new to and from fractions
			- &rarr; computes new local conductance objective value of maximal fraction between from and to via `compute_conductance_objective(..)` of `ConductanceGlobalAttributedGains`  (at least one part isn't empty &rArr; max. isn't $0 / 0$)
			- sets `tmp_scores[to] = local_conductance_after - local_conductance_before`
	- `HyperedgeWeight gain(to_score, isolated_block_gain)`: returns `to_score`

### Part 2.4 Guide: Extending the Library Interface with a Custom Objective Function

#### C Interface

- ```include/mtkahypartypes.h```: Add a enum type to ```mt_kahypar_objective_t``` representing your new objective function:
	+ \+ `CONDUCTANCE_LOCAL`, `CONDUCTANCE_GLOBAL` [`CLUSTER` is also there]
- ```lib/mtkahypar.cpp```: [`CLUSTER` is also there] Create a mapping between the enum types ```mt_kahypar_objective_t``` and ```Objective``` in 
	- ```mt_kahypar_set_context_parameter(...)``` : `"conductance_local"`, `"conductance_global"`
	- ```mt_kahypar_set_partitioning_parameters(...)```: `CONDUCTANCE_LOCAL`, `CONDUCTANCE_GLOBAL`
- ```include/mtkahypar.h```: Add a function that takes a ```mt_kahypar_partitioned_hypergraph_t``` and computes your objective function (similar to ```mt_kahypar_cut(...)``` and ```mt_kahypar_km1```):
	+ \+ `mt_kahypar_conductance_local`, `mt_kahypar_conductance_global` (only declarations)
	+ [my] &rArr; in `mtkahypar.cpp`:
		- \+ implementations of `mt_kahypar_conductance_local`, `mt_kahypar_conductance_global` [analog. to `mt_kahypar_cut(...)`]: both call `lib::conductance_..(phg)` with `Throwing=true`
		- \+ in `mt_kahypar_get_objective(context)` add cases `Objective::conductance_local`, `Objective::conductance_global`
	+ [my] &rarr; `lib_generic_impl.h`: 
		+ \+ `conductance_global(p)`, `conductance_local(p)`: call `metrics::quality` with their `Objective` (analog. to `soed<Throwing>(phg)`)

#### Python Interface

The Python interface is defined in ```python/module.cpp```. 

You only have to add a mapping between a string representation of your new objective function and our ```Objective``` enum type in the enum type section of the file:
- ```python/module.cpp```:
	- Objectives
	```cpp
	using mt_kahypar::Objective;
	py::enum_<Objective>(m, "Objective", py::module_local())
		.value("CUT", Objective::cut)
		.value("KM1", Objective::km1)
		.value("SOED", Objective::soed)
		.value("CONDUCTANCE_LOCAL", Objective::conductance_local)
		.value("CONDUCTANCE_GLOBAL", Objective::conductance_global);
	```
	- preset cluster **!!!**:
	```cpp
  using mt_kahypar::PresetType;
  py::enum_<PresetType>(m, "PresetType", py::module_local())
    .value("DETERMINISTIC", PresetType::deterministic)
    .value("LARGE_K", PresetType::large_k)
    .value("DEFAULT", PresetType::default_preset)
    .value("QUALITY", PresetType::quality)
    .value("HIGHEST_QUALITY", PresetType::highest_quality)
    .value("CLUSTER", PresetType::cluster);
	```

 Afterwards, add function to the ```PartitionedGraph```, ```PartitionedHypergraph``` and ```SparsePartitionedHypergraph``` class that computes the value of your objective function:
 `python/module.cpp`:
 -  PartitionedHypergraph:
```cpp
  // ####################### Partitioned Hypergraph #######################

  phg_class
  	...
    .def("soed", &lib::soed<true>,
      "Computes the sum-of-external-degree metric of the partition")
    .def("conductance_local", &lib::conductance_local<true>,
      "Computes the local conductance metric of the partition")
    .def("conductance_global", &lib::conductance_global<true>,
      "Computes the global conductance metric of the partition")
```