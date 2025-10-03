import os
import matplotlib.pyplot as plt
import numpy as np

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

def plot_data(data_for_letter, data_for_letter_compare, dir_name, scenario_letter, version, version_compare):
    metrics = list(data_for_letter.keys())
    num_metrics = len(metrics)

    plt.figure(figsize=(15, 5 * ((num_metrics + 2) // 3)))
    compare = (data_for_letter_compare is not None)
    if compare:
        plt.suptitle(f'{dir_name}/{scenario_letter} - Version {version} vs {version_compare}', fontsize=16)
    else:
        plt.suptitle(f'{dir_name}/{scenario_letter} - Version {version}', fontsize=16)
    for i, metric in enumerate(metrics):
        data = data_for_letter[metric]
        plt.subplot((num_metrics + 2) // 3, 3, i + 1)
        x_values = []
        y_values = []
        colors = []
        blue = (0.2, 0.4, 0.8, 0.6)
        orange = (1.0, 0.5, 0.0, 0.6)
        black = (0.0, 0.0, 0.0, 0.6)
        
        for number, values in data.items():
            x_values.extend([number] * len(values))
            y_values.extend(values)
            colors.extend([blue] * len(values))
        if compare and metric in data_for_letter_compare:
            data_compare = data_for_letter_compare[metric]
            for number, values in data_compare.items():
                x_values.extend([number] * len(values))
                y_values.extend(values)
                colors.extend([orange] * len(values))
        # display the scattered points in the subplot
        plt.scatter(x_values, y_values, alpha=0.7, c=colors)
        plt.xlabel('Scenario index')
        plt.ylabel(metric)
        plt.title(f'{metric}')

    # add a legende for blue and orange colors
    if compare:
        plt.figlegend([plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=blue, markersize=10, label=f'Version {version}'),
                        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=orange, markersize=10, label=f'Version {version_compare}')],
                       loc='upper right', fontsize=12)
    # Save the plot to a file
    plot_filename = f'plots/{version}/{dir_name}_{scenario_letter}.png'
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
    parser.add_argument("--version", type=str, default="default", help="Version used for the experimental results")
    parser.add_argument("--compare_with_version", type=str, default="none", help="Version to compare with (e.g., AON_HMLL_clique)")
    return parser.parse_args()

# Main function to process all directories and create the plots
def main():
    args = parse_args()
    root_dir = args.data_root_dir
    version = args.version
    compare_with_version = args.compare_with_version
    os.makedirs('plots', exist_ok=True)
    os.makedirs(f"plots/{version}", exist_ok=True)
    
    # Iterate over the directories (DCHSBM, hABCD, HyperSBM)
    for dir_name in os.listdir(root_dir):
        dir_path = os.path.join(root_dir, dir_name)
        
        if os.path.isdir(dir_path):  # Only process directories
            data = process_mode_directory(root_dir, dir_name, version)
            data_compare = dict()
            if compare_with_version != "none":
                data_compare = process_mode_directory(root_dir, dir_name, compare_with_version)

        
        for scenario_letter, data in data.items():
            if scenario_letter in data_compare:
                plot_data(data, data_compare[scenario_letter],
                            dir_name, scenario_letter, version, compare_with_version)
            else:
                plot_data(data, None,
                            dir_name, scenario_letter, version, compare_with_version)

    print("All plots have been generated!")

# Run the main function
if __name__ == "__main__":
    main()
