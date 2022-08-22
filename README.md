# Cho2022_BurstDetection

SungJun Cho & Jee Hyun Choi

Please email Jee Hyun Choi at jeechoi@kist.re.kr or SungJun Cho at scho.sungjun@gmail.com with any questions or concerns.

---

## Getting Started

This repository contains all the scripts necessary to reproduce the analysis and figures shown in our study. To install, simply download this repository folder, specify paths to the folder location, and run the scripts.

The repository is divided into three main categories:

1. `data`: This directory contains the experimental data and simulation results used in this study.
2. `utils`: This directory contains all the functions that are necessary for the analysis and visualization.
3. **Main scripts**: Each `Figure*` or `Table*` directory contains one or more scripts that can be used to reproduce corresponding figures or tables.

## Detailed Description

### Data
* `experimental_data/ESCAPE`
    * `annot`: Contains the `.xlsx` file which includes the starting and ending time points of the bursts annotated by the human expert
    * `lfps`: Contains the `.mat` file that stores LFP signals of the robot-based escape experiment
    * `sample_vidoes`: Contains the sample video recordings of the robot-based escape experiment
* `simulation_data`
    * `HM_*.mat`: The heatmaps of different metrics and algorithms
        * `beta` and `gamma` in the file name indicate that their heatmaps recapitulate algorithmic performances of detecting the beta (23-27 Hz) and gamma (35-45 Hz) bursts, respectively.
    * `randseed_*.mat`: Stores random seeds that were used to simulate synthetic signals
        * `compute_heatmap.m` ensures reproducibility by calling this data and simulating the same signals.
    * `DC_*.mat`: Stores the heatmaps of detection confidence scores
        * These heatmaps were manually computed using the F1-score and temporal concurrence heatmaps in `HM_*.mat`.

### Utility
* `utils` include multiple sub-directories divided by their usage. Every function in the scripts provides a description about its inputs and outputs. Please refer to each script for detail.

### Main Scripts
* Every `.m` script in these directories starts by configuring library paths, which allows an access to `utils` and `data`.
    * **NOTE**: Paths in the scripts are currently set to the ones we used. To run without errors, set them to  the locations where you downloaded the files.
        ```
        util_path = genpath('PATH_TO_UTILS')
        data_path = genpath('PATH_TO_DATA')
        addpath(util_path)
        addpath(data_path)
        ```
* For `.ipynb` and `.py` scripts, you can similarly change the path by setting `input_path`, `file_path`, or `save_path` to your desired location.

### Others
* `compute_heatmap.m`: This script is an important file that should be run prior to executing all the main scripts (except for those in `Figure1` and `Table1`). It outputs simulation results, which are stored as a `struct` format that contains heatmaps for different metrics and algorithms. The heatmaps used in this study are already provided in the `simulation_data` directory, so you do not have to run this script unless you have a specific range of frequency band in which you want to construct the heatmaps.

## Applying the Algorithm Selection Rule
1. To apply the algorithm selecion rule, refer to `fig5C_compute_ecdf.m` and `fig5D_5_2_auc_statistics.m`. These scripts show how you can identify an algorithm optimal for a specific dataset, given the simulation results.
2. As mentioned above, if your frequency of interest differs from the beta and gamma range used in this study, you can compute the heatmaps that suit your interest before applying the selection rule.

## Requirements
The analyses and visualizations in this paper had following dependencies:

```
MATLAB 2020a (or later)
python==3.7.4
seaborn==0.11.1
scipy==1.4.1
numpy==1.17.2
pandas==0.25.1
```
