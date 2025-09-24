/*******************************************************************************
 * MIT License
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

#include "mt-kahypar/partition/initial_partitioning/aon_hypermodularity_bayesian_initial_partitioner.h"

#include "mt-kahypar/definitions.h"
 
 
namespace mt_kahypar {

template<typename TypeTraits>
void AONHypermodularityBayesianInitialPartitioner<TypeTraits>::partitionImpl() {
  /// [debug] std::cout << "AONHypermodularityKernelInitialPartitioner::partitionImpl()" << std::endl;
  _aon_ip.partitionImpl(
    std::max<HypernodeID>(_context.partition.large_hyperedge_size_threshold / 10, 100) /* edgeSizeThreshold */,   
    1e2 /* maxNumIter */,
    1e-8 /* eps */,
    _context.partition.initial_num_nodes / 100 /* clusterPenalty */,
    true /* randomize */,
    true /* useOriginalEdgeSizes */,
    InitialPartitioningAlgorithm::aon_hypermodularity_bayesian /* ip_name */
  );
}

INSTANTIATE_CLASS_WITH_TYPE_TRAITS(AONHypermodularityBayesianInitialPartitioner)
 
} // namespace mt_kahypar