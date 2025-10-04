import os
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# global color vector.
colors = ["blue", "orange", "violet", "purple", "yellow",
           "cyan", "brown", "pink", "gray", "black", "red"]
# green is reserved for ground truth


def read_scenario_data(file_path, letter, number, data_dict):
    ignore_metrics = ["Version", "Instances"]
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
        # Loop through each line in the file and parse it
        for line in lines:
            if '=' in line:
                metric_name, values_str = line.split('=', 1)
                metric_name = metric_name.strip()
                if metric_name in ignore_metrics:
                    continue
                # Strip brackets and split the numbers into a list
                values_str = values_str.strip().strip('[]')
                values = values_str.split(',')
                values = [float(v.strip()) for v in values if v.lower() != "nan"]
                if metric_name not in data_dict[letter]:
                    data_dict[letter][metric_name] = dict()
                data_dict[letter][metric_name][number] = values

def get_points_and_means(data, metric):
    x_values = []
    y_values = []
    colors = []
    blue = (0.2, 0.4, 0.8)
    orange = (1.0, 0.5, 0.0)
    mean_points = []
    
    for number, values in data.items():
        x_values.extend([number] * len(values))
        y_values.extend(values)
        # add mean point
        mean = np.mean(values)
        mean_points.append((number, mean))
    return x_values, y_values, mean_points



def plot_data(version_to_data, versions, dir_name, scenario_letter, fout_name, colors = colors):
    gt_for_metric = {
        "K": "GT_K",
        "Modularity": "GT_Mod",
        "Conductance": "GT_Conductance",
    }
    skip_metrics = gt_for_metric.values()
    
    metrics_set = set()
    for version in versions:
        metrics_set.update(version_to_data[version][scenario_letter].keys())
    metrics = [m for m in sorted(metrics_set) if m not in skip_metrics]
    num_metrics = len(metrics)
    has_gt = False

    plt.figure(figsize=(15, 5 * ((num_metrics + 2) // 3)))
    plt.suptitle(f'{dir_name}/{scenario_letter} - Versions: {", ".join(versions)}', fontsize=16)
    for j, version in enumerate(versions):
        letter_to_data = version_to_data[version]
        print(f"Processing version: {version} number {j} with data for letters: {list(letter_to_data.keys())}")
        version_color = colors[j]
        metric_to_data = letter_to_data.get(scenario_letter, {})
        for i, metric in enumerate(metrics):
            if metric in skip_metrics:
                continue
            if metric not in metric_to_data:
                continue
            data = metric_to_data[metric]
            plt.subplot((num_metrics + 2) // 3, 3, i + 1)

            
            # plot ground truth if available
            if metric in gt_for_metric.keys():
                print(f"Found ground truth for metric {metric}")
                gt_metric = gt_for_metric[metric]
                if gt_metric in metric_to_data:
                    gt_data = metric_to_data[gt_metric]
                    # plot the mean ground truth line
                    _, _, gt_mean_points = get_points_and_means(gt_data, gt_metric)
                    if gt_mean_points:
                        has_gt = True
                        gt_mean_points.sort(key=lambda x: x[0])
                        mean_x, mean_y = zip(*gt_mean_points)
                        plt.plot(mean_x, mean_y, color='green', linewidth=2, label=f'Mean {gt_metric}', linestyle='--', marker='o', alpha=0.3)

            # plot version data
            x_values, y_values, mean_points = get_points_and_means(data, metric)
            plt.scatter(x_values, y_values, c=version_color, alpha=0.4, edgecolors='gray')
            if mean_points:
                mean_points.sort(key=lambda x: x[0])
                mean_x, mean_y = zip(*mean_points)
                plt.plot(mean_x, mean_y, color=version_color, linewidth=2, label=f'Mean {version}', linestyle='--', marker='_', alpha=0.5)
            

            plt.xlabel('Scenario index')
            plt.ylabel(metric)
            plt.title(f'{metric}')

        # add a legende
        patch_handles = []
        for version in versions:
            version_color = colors[versions.index(version)]
            patch = mpatches.Patch(color=version_color, label=version, linewidth=2)
            patch_handles.append(patch)
        if has_gt:
            green_patch = mpatches.Patch(color='green', label=f'Ground Truth (Mean)', linewidth=2)
            patch_handles.append(green_patch)

        plt.figlegend(handles=patch_handles, loc='upper right', fontsize=12)
        # Save the plot to a file
        # plot_filename = f'plots/{dir_name}_{scenario_letter}_{",".join(versions)}.png'
        plot_filename = f'plots/{fout_name}/{dir_name}_{scenario_letter}.png'
        plt.savefig(plot_filename)
        print(f"Plot saved to {plot_filename}")

# Function to process a mode directory
def process_mode_directory(root_dir, dir_name, version):
    data_dict = dict()
    for subdir in os.listdir(root_dir + '/' + dir_name):
        subdir_path = os.path.join(root_dir, dir_name, subdir)
        
        if os.path.isdir(subdir_path):  # Only process directories
            # Check format scenXk
            if subdir[0:4] != 'scen' or not subdir[4].isalpha() or not subdir[5:].isdigit():
                continue
            scenario_letter, scenario_number = subdir[4], int(subdir[5:])
            data_file = f"{subdir_path}/{version}.results"
            if not os.path.isfile(data_file):
                print(f"Warning: {data_file} does not exist. Skipping.")
                continue

            # save the data to data_dict
            if scenario_letter not in data_dict:
                data_dict[scenario_letter] = dict()
            read_scenario_data(data_file, scenario_letter, scenario_number, data_dict)
    ######################
    return data_dict


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Plot clustering metrics from experimental results.")
    parser.add_argument("--data_root_dir", type=str, help="Root directory containing the mode subdirectories")
    parser.add_argument("--versions", type=str, default="default", help="Comma-separated versions to be used for plotting")
    parser.add_argument("--output", type=str, default="plot", help="Base name for output plot files")
    return parser.parse_args()

# Main function to process all directories and create the plots
def main():
    args = parse_args()
    root_dir = args.data_root_dir
    versions = args.versions.split(",")
    fout_name = args.output
    os.makedirs('plots', exist_ok=True)
    if fout_name != "":
        os.makedirs(f'plots/{fout_name}', exist_ok=True)

    if len(versions) > len(colors):
        print(f"Error: Number of versions ({len(versions)}) exceeds number of available colors ({len(colors)}).)")
        print("Please add more colors to the colors list in the script plot_all_data.py.")
        return

    # Iterate over the directories (DCHSBM, hABCD, HyperSBM)
    for dir_name in os.listdir(root_dir):
        dir_path = os.path.join(root_dir, dir_name)

        letters = set()
        
        if os.path.isdir(dir_path):  # Only process directories
            version_to_data = dict()
            for version in versions:
                version_to_data[version] = process_mode_directory(root_dir, dir_name, version)
                letters.update(version_to_data[version].keys())

        for scenario_letter in letters:
            plot_data(version_to_data, versions,
                        dir_name, scenario_letter, fout_name)

    print("All plots have been generated!")

# Run the main function
if __name__ == "__main__":
    main()
