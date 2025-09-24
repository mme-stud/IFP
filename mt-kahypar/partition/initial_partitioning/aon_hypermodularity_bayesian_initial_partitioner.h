/*******************************************************************************
 * MIT License
 * 
 * This file is part of Mt-KaHyPar.
 *
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

#pragma once

#include "mt-kahypar/partition/initial_partitioning/i_initial_partitioner.h"
#include "mt-kahypar/partition/initial_partitioning/initial_partitioning_data_container.h"
#include "mt-kahypar/parallel/stl/scalable_vector.h"
#include "mt-kahypar/datastructures/array.h"

#include "mt-kahypar/partition/initial_partitioning/aon_hypermodularity_initial_partitioner.h"

namespace mt_kahypar {

template<typename TypeTraits>
class AONHypermodularityBayesianInitialPartitioner : public IInitialPartitioner {

  static constexpr bool debug = false;

  using PartitionedHypergraph = typename TypeTraits::PartitionedHypergraph;
  using UnderlyingHypergraph = typename PartitionedHypergraph::UnderlyingHypergraph;

 public:
  AONHypermodularityBayesianInitialPartitioner(const InitialPartitioningAlgorithm,
                           ip_data_container_t* ip_data,
                           const Context& context,
                           const int seed, const int tag) :
    _aon_ip(InitialPartitioningAlgorithm::aon_hypermodularity_bayesian, ip_data, context, seed, tag),
    _context(context)
    { }

 private:
  
  void partitionImpl() final;

  AONHypermodularityInitialPartitioner<TypeTraits> _aon_ip;
  const Context& _context;
};

} // namespace mt_kahypar
