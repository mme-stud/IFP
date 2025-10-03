#include <cassert>
#include <memory>
#include <vector>
#include <thread>

#include <string>
#include <iostream>
#include <fstream>

#include <tuple>

#include <mtkahypar.h>

std::string survey_root_dir = "survey_benchmark";

std::tuple<std::string, std::string, std::string, int, std::vector<int>>
parse_commandline_args(int argc, char* argv[]) {
  std::string version = "default"; // code version / configuration
  std::string mode = "none"; // DCHSBM, hABCD, HyperSBM
  std::string scenario = "none"; // scenA1... 
  std::string graphs = "none"; // comma-separated list of graph ids, e.g. "1,2,3"
  std::vector<int> instances; // parsed from graphs

  int num_threads = -1;
  if (argc != 1 + 10 && argc != 1 + 8) {
    std::cout << "Usage: -v <version> -m HyperSBM -s scenA1 -t <num_threads> [-instances \"1,22\"]" << std::endl;
    std::exit(1);
  }
  for (int i = 1; i < argc; i+=2) {
    std::string arg = argv[i];
    if (arg == "-v") {
      version = argv[i+1];
    } else if (arg == "-m") {
      mode = argv[i+1];
    } else if (arg == "-s") {
      scenario = argv[i+1];
    } else if (arg == "-t") {
      num_threads = std::stoi(argv[i+1]);
    } else if (arg == "-instances") {
      // comma-separated list of graph ids
      graphs = argv[i+1];
    } else {
      std::cout << "Unknown command line argument: " << arg << std::endl;
      std::exit(1);
    }
  }
  if (graphs == "none") {
    std::cout << "Graphs not specified, using default: all" << std::endl;
    for (int i = 1; i <= 25; ++i) {
      graphs = graphs + std::to_string(i) + ","; 
    }
    graphs.pop_back(); // remove last comma
  }
  // parse graphs
  size_t pos = 0;
  while ((pos = graphs.find(',')) != std::string::npos) {
    instances.push_back(std::stoi(graphs.substr(0, pos)));
    graphs.erase(0, pos + 1);
  }
  if (!graphs.empty()) {
    instances.push_back(std::stoi(graphs));
  }
  // output parsed args
  std::cout << "Version: " << version << std::endl;
  std::cout << "Mode: " << mode << std::endl;
  std::cout << "Scenario: " << scenario << std::endl;
  std::cout << "Number of Threads: " << num_threads << std::endl;
  std::cout << "Instances: ";
  for (const auto& inst : instances) {
    std::cout << inst << " ";
  }
  std::cout << std::endl;

  assert(num_threads > 0);
  assert(num_threads <= std::thread::hardware_concurrency());
  return {version, mode, scenario, num_threads, instances};
}


std::tuple<double, double> 
run_partitioning(const std::string& hgr_path, const std::string& hgr_name, mt_kahypar_error_t& error, int num_threads, std::string version) {
  std::cout << "Hypergraph: " << hgr_path << "/" << hgr_name << std::endl;
  mt_kahypar_initialize(num_threads, true /* activate interleaved NUMA allocation policy */);
  // Setup partitioning context
  mt_kahypar_context_t* context = mt_kahypar_context_from_preset(CLUSTER);
  mt_kahypar_set_partitioning_parameters(context,
    2 /* dummy number of blocks */, 0.03 /* dummy imbalance parameter */,
    CONDUCTANCE_LOCAL /* objective function */);
  mt_kahypar_set_seed(42 /* seed */);

  // Disable logging
  mt_kahypar_status_t status =
    mt_kahypar_set_context_parameter(context, VERBOSE, "0", &error);
  assert(status == SUCCESS);
  
  std::cout << "Input: start" << std::endl;
  // time input
  auto startTime = std::chrono::high_resolution_clock::now();
  // Load Hypergraph for CLUSTER preset
  const std::string hgr_file = hgr_path + ("/") + hgr_name;
  mt_kahypar_hypergraph_t hypergraph =
    mt_kahypar_read_hypergraph_from_file(hgr_file.c_str(),
      context, HMETIS /* file format */, &error);
  auto endTime = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration<double>(endTime - startTime);
  std::cout << "Time:  " << duration.count() << " seconds.\n";
  if (hypergraph.hypergraph == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }
  std::cout << "Input: end" << std::endl;

  // Partition Hypergraph
  std::cout << "Partition: start" << std::endl;
  // time partition
  startTime = std::chrono::high_resolution_clock::now();
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_partition(hypergraph, context, &error);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration<double>(endTime - startTime);
  std::cout << "Time:  " << duration.count() << " seconds.\n";
  std::cout << "Partition: end" << std::endl;
  if (partitioned_hg.partitioned_hg == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }
 
  // Output some info about the partition
  { 
    // Extract Partition
    auto partition = std::make_unique<mt_kahypar_partition_id_t[]>(
      mt_kahypar_num_hypernodes(hypergraph));
    mt_kahypar_get_partition(partitioned_hg, partition.get());
    // k - not from context [context is copied internaly to mt-kahypar and not updated]
    int k = mt_kahypar_num_blocks(partitioned_hg);
    // Extract Block Weights
    auto block_weights = std::make_unique<mt_kahypar_hypernode_weight_t[]>(k);
    mt_kahypar_get_block_weights(partitioned_hg, block_weights.get());

    // Get the number of nonempty blocks
    mt_kahypar_partition_id_t num_nonempty_blocks = 0;
    for (int i = 0; i < k; ++i) {
      if (block_weights[i] > 0) {
        ++num_nonempty_blocks;
      }
    }
    std::cout << "Number of clusters: " << num_nonempty_blocks << std::endl;
    std::cout << "Number of blocks: " << k << std::endl;
    std::cout << "Computed conductance: " << mt_kahypar_conductance_local(partitioned_hg) << std::endl;
    std::cout << "Computed AON hypermodularity: " << mt_kahypar_aon_hypermodularity(partitioned_hg) << std::endl;
  }

  // Write the partition to a file
  const std::string partition_file = hgr_path + "/" + version + "/" + hgr_name + ".part";
  mt_kahypar_write_partition_to_file(
        partitioned_hg, 
        partition_file.c_str(),
        &error);
  std::cout << "Wrote partition to file: " << partition_file << std::endl;

  double time = duration.count();
  double conductance = mt_kahypar_conductance_local(partitioned_hg);

  mt_kahypar_free_context(context);
  mt_kahypar_free_hypergraph(hypergraph);
  mt_kahypar_free_partitioned_hypergraph(partitioned_hg);

  return {time, conductance};
}

// args: mode, scenario, num_threads, [num_graphs]
int main(int argc, char* argv[]) {
  mt_kahypar_error_t error{}; // to capture error messages from mt-kahypar functions
  
  // parse command line args
  auto parsed_args = parse_commandline_args(argc, argv);
  std::string version = std::get<0>(parsed_args);
  std::string mode = std::get<1>(parsed_args);
  std::string scenario = std::get<2>(parsed_args);
  int num_threads = std::get<3>(parsed_args);
  std::vector<int> instances = std::get<4>(parsed_args);
  int num_graphs = instances.size();

  std::vector<double> times(num_graphs + 1, 0);
  std::vector<double> conductances(num_graphs + 1, 0);
  
  // Loop through the graphs in the scenario
  std::string graph_path = survey_root_dir + "/" + mode + "/" + scenario;
  for (int i = 1; i <= num_graphs; ++i) {
    int graph_id = instances[i - 1];
    std::string graph_name = std::string("rep") + (std::to_string(graph_id)) + ("_he.hgr");
    // write the partition to the file mt_kahypar_<graph_filename>_<num_threads>threads.part
    auto time_conductance = run_partitioning(graph_path, graph_name, error, num_threads, version);
    times[i] = std::get<0>(time_conductance);
    conductances[i] = std::get<1>(time_conductance);
  }
  // Write all times and conductances to a file
  std::string results_filename = 
        survey_root_dir + ("/") + mode + ("/") 
      + scenario + ("/") + version + (".results.t_c");
  std::ofstream results_file(results_filename);

  results_file << "Version = " << version << "\n";
  
  results_file << "Instances = [";
  for (int i = 1; i <= num_graphs; ++i) {
    results_file << instances[i - 1];
    if (i < num_graphs) {
      results_file << ", ";
    }
  }
  results_file << "]\n";
  
  results_file << "Time (sec.) = [";
  for (int i = 1; i <= num_graphs; ++i) {
    results_file << times[i];
    if (i < num_graphs) {
      results_file << ", ";
    }
  }
  results_file << "]\n";
  results_file << "Conductance = [";
  for (int i = 1; i <= num_graphs; ++i) {
    results_file << conductances[i];
    if (i < num_graphs) {
      results_file << ", ";
    }
  }
  results_file << "]\n";
  results_file.close();
  std::cout << "Wrote times and conductances to file: " << results_filename << std::endl;
  return 0;
}