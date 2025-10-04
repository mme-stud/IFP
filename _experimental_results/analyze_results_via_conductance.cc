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
std::string preset_dir = "config/";

std::tuple<std::string, std::string, std::string, std::string, std::vector<int>>
parse_commandline_args(int argc, char* argv[]) {
  std::string version = "default"; // code version / configuration
  std::string preset = "cluster"; // mt-kahypar preset
  std::string mode = "none"; // DCHSBM, hABCD, HyperSBM
  std::string scenario = "none"; // scenA1... 
  std::string graphs = "none"; // comma-separated list of graph ids, e.g. "1,2,3"
  std::vector<int> instances; // parsed from graphs

  int num_threads = -1;
  if (argc != 1 + 10 && argc != 1 + 8) {
    std::cout << "[Analyze conductance] Usage: -p <preset> -v <version> -m HyperSBM -s scenA1  [-instances \"1,22\"]" << std::endl;
    std::exit(1);
  }
  for (int i = 1; i < argc; i+=2) {
    std::string arg = argv[i];
    if (arg == "-p") {
      preset = argv[i+1];
    } else if (arg == "-v") {
      version = argv[i+1];
    } else if (arg == "-m") {
      mode = argv[i+1];
    } else if (arg == "-s") {
      scenario = argv[i+1];
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
  std::cout << "Parsed command line arguments:" << std::endl;
  std::cout << "Preset: " << preset << std::endl;
  std::cout << "Version: " << version << std::endl;
  std::cout << "Mode: " << mode << std::endl;
  std::cout << "Scenario: " << scenario << std::endl;
  std::cout << "Instances: ";
  for (const auto& inst : instances) {
    std::cout << inst << " ";
  }
  std::cout << std::endl;

  return {preset, version, mode, scenario, instances};
}


double 
compute_conductance(
        const std::string& hgr_path, const std::string& hgr_name, 
        const std::string& partition_file,  
        mt_kahypar_error_t& error, std::string preset) {
  std::cout << "Hypergraph: " << hgr_path << "/" << hgr_name << std::endl;
  mt_kahypar_initialize(1, true /* activate interleaved NUMA allocation policy */);
  // Setup partitioning context from preset
  mt_kahypar_context_t* context = mt_kahypar_context_from_file((preset_dir + preset + "_preset.ini").c_str(), &error);
  if (context == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }
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
  std::cout << "Partition construction: start" << std::endl;
  // time partition
  startTime = std::chrono::high_resolution_clock::now();
  // read number of blocks from the partition file
  std::ifstream pf(partition_file);
  if (!pf) {
    std::cout << "Partition file " << partition_file << " not found." << std::endl; std::exit(1);
  }
  int num_blocks = 0;
  std::string line;
  while (std::getline(pf, line)) {
    int block = std::stoi(line);
    if (block + 1 > num_blocks) {
      num_blocks = block + 1;
    }
  }
  // read partition from file
  mt_kahypar_partitioned_hypergraph_t partitioned_hg =
    mt_kahypar_read_partition_from_file(hypergraph,
                                        context,
                                        num_blocks,
                                        partition_file.c_str(),
                                        &error);
  endTime = std::chrono::high_resolution_clock::now();
  duration = std::chrono::duration<double>(endTime - startTime);
  std::cout << "Time:  " << duration.count() << " seconds.\n";
  std::cout << "Partition: end" << std::endl;
  if (partitioned_hg.partitioned_hg == nullptr) {
    std::cout << error.msg << std::endl; std::exit(1);
  }
 
  double conductance = mt_kahypar_conductance_local(partitioned_hg);
  std::cout << "Computed conductance: " << conductance << std::endl;
  mt_kahypar_free_context(context);
  mt_kahypar_free_hypergraph(hypergraph);
  mt_kahypar_free_partitioned_hypergraph(partitioned_hg);
  return conductance;
}

// args: mode, scenario, num_threads, [num_graphs]
int main(int argc, char* argv[]) {
  mt_kahypar_error_t error{}; // to capture error messages from mt-kahypar functions
  
  // parse command line args
  auto parsed_args = parse_commandline_args(argc, argv);
  std::string preset = std::get<0>(parsed_args);
  std::string version = std::get<1>(parsed_args);
  std::string mode = std::get<2>(parsed_args);
  std::string scenario = std::get<3>(parsed_args);
  std::vector<int> instances = std::get<4>(parsed_args);
  int num_graphs = instances.size();

  std::vector<double> conductances_hat(num_graphs + 1, 0);
  std::vector<double> conductances_true(num_graphs + 1, 0);

  
  // Loop through the graphs in the scenario
  std::string graph_path = survey_root_dir + "/" + mode + "/" + scenario;
  for (int i = 1; i <= num_graphs; ++i) {
    int graph_id = instances[i - 1];
    std::string graph_name = std::string("rep") + (std::to_string(graph_id)) + ("_he.hgr");
    std::string partition_file = graph_path + ("/") + version + ("/") + graph_name + ".part";
    std::string gt_file = graph_path + ("/rep") + (std::to_string(graph_id)) + ("_assign.txt");
    // write the partition to the file mt_kahypar_<graph_filename>_<num_threads>threads.part
    conductances_hat[i] = compute_conductance(graph_path, graph_name, partition_file, error, preset);
    conductances_true[i] = compute_conductance(graph_path, graph_name, gt_file, error, preset);
  }
  // Write all times and conductances to a file
  std::string results_filename = 
        survey_root_dir + ("/") + mode + ("/") 
      + scenario + ("/") + version + (".results");
  // APPEND to existing file
  std::ofstream results_file(results_filename, std::ios_base::app);

  results_file << "Conductance = [";
  for (int i = 1; i <= num_graphs; ++i) {
    results_file << conductances_hat[i];
    if (i < num_graphs) {
      results_file << ", ";
    }
  }
  results_file << "]\n";

  results_file << "GT_Conductance = [";
  for (int i = 1; i <= num_graphs; ++i) {
    results_file << conductances_true[i];
    if (i < num_graphs) {
      results_file << ", ";
    }
  }
  results_file << "]\n";

  results_file.close();
  std::cout << "Wrote times and conductances to file: " << results_filename << std::endl;
  return 0;
}