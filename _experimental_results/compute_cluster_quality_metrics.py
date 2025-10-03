import numpy as np
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score, normalized_mutual_info_score, homogeneity_score, completeness_score, v_measure_score
from collections import Counter
import sys

def read_clusters_from_file(file_path):
    """
    Reads cluster IDs from a file where each line corresponds to a node's cluster ID.
    """
    try:
        with open(file_path, "r") as file:
            return [int(line.strip()) for line in file]
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
        sys.exit(1)

def mutual_information(Z, Z_hat, normalized=False):
    """
    Compute the mutual information between two clusterings (Z and Z_hat).
    Optionally normalize to calculate NMI.
    """
    n = len(Z)

    # Joint probabilities
    p_XY = Counter((Z[i], Z_hat[i]) for i in range(n))
    for key in p_XY:
        p_XY[key] /= n

    # Marginal probabilities
    p_X = Counter(Z)  # Marginal for Z
    for key in p_X:
        p_X[key] /= n

    p_Y = Counter(Z_hat)  # Marginal for Z_hat
    for key in p_Y:
        p_Y[key] /= n

    # Mutual Information
    I = 0
    for (x, y), p_xy in p_XY.items():
        try:
            I += p_xy * np.log(p_xy / (p_X[x] * p_Y[y]))
        except (ValueError, KeyError):
            pass  # Handle log(0) or missing keys gracefully

    if normalized:
        # Compute entropies
        H_X = -sum(p * np.log(p) for p in p_X.values())
        H_Y = -sum(p * np.log(p) for p in p_Y.values())
        return 2 * I / (H_X + H_Y) if (H_X + H_Y) > 0 else 0

    return I

def compute_purity(predicted_clusters, ground_truth_clusters):
    """
    Compute purity for clustering.
    """
    total_samples = len(ground_truth_clusters)
    contingency_matrix = Counter(zip(predicted_clusters, ground_truth_clusters))

    max_matching = 0
    for cluster in set(predicted_clusters):
        max_matching += max(
            contingency_matrix.get((cluster, gt), 0)
            for gt in set(ground_truth_clusters)
        )
    return max_matching / total_samples

def compute_f1_score(predicted_clusters, ground_truth_clusters):
    """
    Compute F1 score for clustering using pairwise precision and recall.
    """
    n = len(predicted_clusters)

    # Count pairs of points that are in the same cluster
    def count_same_cluster_pairs(clusters):
        cluster_count = Counter(clusters)
        return sum(v * (v - 1) // 2 for v in cluster_count.values())

    # Calculate true positives, precision, and recall
    true_positive = sum(
        1 for i in range(n) for j in range(i + 1, n)
        if (predicted_clusters[i] == predicted_clusters[j] and ground_truth_clusters[i] == ground_truth_clusters[j])
    )
    predicted_positive = count_same_cluster_pairs(predicted_clusters)
    actual_positive = count_same_cluster_pairs(ground_truth_clusters)

    precision = true_positive / predicted_positive if predicted_positive > 0 else 0
    recall = true_positive / actual_positive if actual_positive > 0 else 0

    # Calculate F1 score
    if precision + recall == 0:
        return 0.0
    return 2 * (precision * recall) / (precision + recall)

def compute_clustering_metrics(predicted_clusters, ground_truth_clusters):
    """
    Compute various clustering quality metrics.

    Returns:
        dict: A dictionary with all computed metrics.
    """
    metrics = {
        "AMI": adjusted_mutual_info_score(ground_truth_clusters, predicted_clusters),
        "ARI": adjusted_rand_score(ground_truth_clusters, predicted_clusters),
        "NMI (sklearn)": normalized_mutual_info_score(ground_truth_clusters, predicted_clusters),
        "NMI (custom)": mutual_information(ground_truth_clusters, predicted_clusters, normalized=True),
        "Mutual Information (custom)": mutual_information(ground_truth_clusters, predicted_clusters),
        "Homogeneity": homogeneity_score(ground_truth_clusters, predicted_clusters),
        "Completeness": completeness_score(ground_truth_clusters, predicted_clusters),
        "V-Measure": v_measure_score(ground_truth_clusters, predicted_clusters),
        "Purity": compute_purity(predicted_clusters, ground_truth_clusters),
        "F1 Score (Pairwise)": compute_f1_score(predicted_clusters, ground_truth_clusters),
    }
    return metrics

def print_cluster_info(clusters, label="Clusters"):
    print(f"\n--- {label} Info ---")
    counter = Counter(clusters)
    print(f"Number of clusters: {len(counter)}")
    print("Cluster sizes:")
    for cid, size in sorted(counter.items()):
        print(f"  Cluster {cid}: {size} nodes")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <predicted_clusters_file> <ground_truth_clusters_file>")
        sys.exit(1)

    predicted_file = sys.argv[1]
    ground_truth_file = sys.argv[2]

    # Read clusters from files
    predicted_clusters = read_clusters_from_file(predicted_file)
    ground_truth_clusters = read_clusters_from_file(ground_truth_file)

    print_cluster_info(predicted_clusters, label="Predicted")
    print_cluster_info(ground_truth_clusters, label="Ground Truth")

    # Compute and print clustering metrics
    scores = compute_clustering_metrics(predicted_clusters, ground_truth_clusters)
    for metric, score in scores.items():
        print(f"{metric}: {score:.4f}")


'''
Clustering Quality Metrics

    Adjusted Mutual Information (AMI): Adjusts for chance when comparing mutual information between two cluster assignments. Higher values indicate better clustering.

    Adjusted Rand Index (ARI): Measures the similarity between two cluster assignments, corrected for chance. A perfect score is 1, and random clustering gets an expected ARI near 0.

    Normalized Mutual Information (NMI): Measures the amount of information shared between two clusterings, normalized to scale from 0 to 1. Unlike AMI, it is not adjusted for chance.

    Purity: Measures how well a clustering matches the ground truth by associating each cluster with the most frequent ground-truth label. Higher purity indicates better clustering.

    F1 Measure: A harmonic mean of precision and recall for clustering. Precision measures the proportion of predicted clusters that are relevant, and recall measures the proportion of relevant clusters that are predicted.

    Homogeneity: A cluster satisfies homogeneity if all its members belong to the same class.

    Completeness: A cluster satisfies completeness if all members of a given class are assigned to the same cluster.

    V-Measure: The harmonic mean of homogeneity and completeness.
'''