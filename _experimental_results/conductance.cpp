#include <cassert>
#include <memory>
#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>


std::string path_to_file;
std::string path_to_parts;
bool mt_kahypar;

void computeConductance(std::string filename_hg, std::string filename_part) {
  std::cout << "Dataset: " << filename_hg << std::endl;
  // open file with edges
  std::string path_hg = path_to_file + filename_hg;
  std::ifstream input_hg(path_hg);
  if (!input_hg) {
    std::cerr << "Error opening file: " << path_hg << std::endl;
    return;
  }
  long long num_edges, num_nodes;
  input_hg >> num_edges >> num_nodes;
  std::cout << "Num edges: " << num_edges << std::endl;
  std::cout << "Num nodes: " << num_nodes << std::endl;

  // open file with partitions
  std::string path_part = path_to_parts + filename_part;
  std::ifstream input_parts(path_part);
  if (!input_parts) {
    std::cerr << "Error opening file: " << path_part << std::endl;
    return;
  }
  // read parts
  std::vector<long long> partID(num_nodes, -1);
  std::cout << "Input parts: start" << std::endl;
  long long num_parts = 0;
  long long node = 0;
  for( std::string line; std::getline( input_parts, line ); ) {
    std::istringstream iss(line);
    while (iss >> partID[node]) {
      num_parts = std::max(num_parts, partID[node] + 1);
      node++;
      if (node > num_nodes) {
        std::cerr << "Error: too many nodes in partition file" << std::endl;
        return;
      }      
    }
  }
  if (node != num_nodes) {
    std::cerr << "Error: not enough nodes in partition file" << std::endl;
    return;
  }
  // close file
  input_parts.close();
  std::cout << "Num parts: " << num_parts << std::endl;
  std::cout << "Input parts: end" << std::endl;

  // read edges
  std::vector<long long> part_volume(num_parts, 0);
  std::vector<long long> cut_weight(num_parts, 0);
  long long total_volume = 0; 
  std::cout << "Input edges: start" << std::endl;
  for( std::string line; std::getline( input_hg, line ); ) {
    std::istringstream iss(line);
    std::vector<bool> cuts(num_parts, false);
    long long num_blocks = 0;
    long long node = 0;
    while (iss >> node) {
      if (mt_kahypar)
        node --;
      if (node >= num_nodes) {
        std::cerr << "Error: too many nodes in edge file" << std::endl;
        return;
      }
      long long part_node = partID[node];
      if (part_node == -1) {
        std::cout << node << " has part " << -1 << std::endl;
        return;
      }
      total_volume++;
      part_volume[part_node]++;
      if ( !cuts[part_node] ) {
        cuts[part_node] = true;
        num_blocks++;
      }
    }
    if ( num_blocks > 1 ) {
      for (long long i = 0; i < num_parts; ++i) {
        if (cuts[i]) {
          cut_weight[i]++;
        }
      }
    }
  }
  // close file
  input_hg.close();

  double conductance = 0.0;
  // std::cout << total_volume << std::endl;
  for (long long i = 0; i < num_parts; ++i) {
    // std::cout << cut_weight[i] << " " << part_volume[i] << std::endl;
    if (part_volume[i] == 0) {
      continue;
    }
    if (part_volume[i] == total_volume) {
      std::cout << "Only one cluster" << std::endl;
      break;
    }
    double c = static_cast<double>(cut_weight[i]) /
                static_cast<double>(std::min(part_volume[i], total_volume - part_volume[i]));
    if (c >= conductance) {
      conductance = c;
    }
  }
  std::cout << "Conductance: " << conductance << std::endl;
  std::cout << "=============================" << std::endl;
}


int main(int argc, char* argv[]) {
  std::cout << "Path to graph file: ";
  std::cin >> path_to_file;
  std::cout << "Graph file: ";
  std::string filename_hg;
  std::cout << "Path to parts: ";
  std::cin >> path_to_parts;
  std::cout << "Parts file: ";
  std::string filename_part;
  std::cout << "Is is mt_kahypar? (0/1): ";
  std::cin >> mt_kahypar;
  std::cout << "=============================" << std::endl;
  computeConductance(filename_hg, filename_part);
}