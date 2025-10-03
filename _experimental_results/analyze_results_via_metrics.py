# use compute_cluster_quality_metrics.py to compute clustering metrics
import os
import sys
import argparse
from compute_cluster_quality_metrics import read_clusters_from_file, print_cluster_info, compute_clustering_metrics

def get_scores(predicted_file, ground_truth_file):    
    # Read clusters from files
    predicted_clusters = read_clusters_from_file(predicted_file)
    ground_truth_clusters = read_clusters_from_file(ground_truth_file)

    print_cluster_info(predicted_clusters, label="Predicted")
    print_cluster_info(ground_truth_clusters, label="Ground Truth")

    # Compute and print clustering metrics
    scores = compute_clustering_metrics(predicted_clusters, ground_truth_clusters)
    return scores

def add_results(results, graph_id, scores):
    for metric, score in scores.items():
        if metric not in results:
            results[metric] = []
        results[metric].append(score)

def parse_args():
    parser = argparse.ArgumentParser(description="Analyze clustering results via various metrics.")
    parser.add_argument("--version", type=str, default="default", help="Version identifier for the experiment")
    parser.add_argument("--data_root_dir", type=str, default = "survey_benchmark", help="Root directory containing the mode isrectoried")
    parser.add_argument("--mode", type=str, help="DCHSBM, hABCD or HyperSBM")
    parser.add_argument("--scenario", type=str, help="e.g. scenA1")
    parser.add_argument("--graphs", type=str, help="Comma-separated list of graph indices to analyze, e.g. 1,2,5")
    return parser.parse_args()


# Usage example:
# python script.py --version <version> --data_root_dir /path/to/data --mode DCHSBM  --scenario scenA1 --graphs 1,2,5 

# parse command line arguments
args = parse_args()
version = args.version
data_root_dir = args.data_root_dir
mode = args.mode
scenario = args.scenario
graphs = args.graphs

output_file = f"{data_root_dir}/{mode}/{scenario}/{version}.results"
results = dict()

for graph_id in graphs.split(","):
    predicted_file = f"{data_root_dir}/{mode}/{scenario}/{version}/rep{graph_id}_he.hgr.part"
    ground_truth_file = f"{data_root_dir}/{mode}/{scenario}/rep{graph_id}_assign.txt"
    print(f"\nAnalyzing graph {graph_id} in mode {mode}, scenario {scenario}")
    scores = get_scores(predicted_file, ground_truth_file)
    add_results(results, graph_id, scores)

# Write results to output file
metrics = ["ARI", "NMI (custom)", "Purity", "F1 Score (Pairwise)"]

# append results, if the file already exists
# if the file does not exist, create it
if not os.path.exists(output_file):
    with open(output_file, "w") as f:
        f.write(f"Version = {version}\n")
        for metric in metrics:
            if metric in results:
                scores = results[metric]
                scores_str = ", ".join(f"{score:.4f}" for score in scores)
                f.write(f"{metric} = [{scores_str}]\n")
else:
    with open(output_file, "a") as f:
        f.write(f"\nVersion = {version}\n")
        for metric in metrics:
            if metric in results:
                scores = results[metric]
                scores_str = ", ".join(f"{score:.4f}" for score in scores)
                f.write(f"{metric} = [{scores_str}]\n")
                print("output_file:", output_file)