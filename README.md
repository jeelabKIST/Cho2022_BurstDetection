# Cho2022_BurstDetection

🔗 A guide towards optimal detection of transient oscillatory bursts with unknown parameters (DOI: [10.1088/1741-2552/acdffd](https://doi.org/10.1088/1741-2552/acdffd))

💡 Please email Jee Hyun Choi at jeechoi@kist.re.kr or SungJun Cho at scho.sungjun@gmail.com with any questions or concerns.

---

## ⚡️ Getting Started

This repository contains all the scripts necessary to reproduce the analysis and figures shown in our study. To install, simply download this repository folder, specify paths to the folder location, and run the scripts.

The repository is divided into three main categories:

1. `data`: This directory contains the experimental data and simulation results used in this study.
2. `utils`: This directory contains all the functions that are necessary for the analysis and visualization.
3. `tutorials`: This directory contains user guidelines for readers who wish to reproduce the results or implement our algorithm selection rule.
4. **Main scripts**: Each `Figure*` or `Table*` directory contains one or more scripts that can be used to reproduce corresponding figures or tables.

## 📄 Detailed Descriptions

### Data
* Experimental Data

   | Directory                                        | Description                                                                               |
   | :----------------------------------------------- | :---------------------------------------------------------------------------------------- |
   | [annot](https://github.com/jeelabKIST/Cho2022_BurstDetection/tree/main/data/experimental_data/ESCAPE/annot)| Contains a `.xlsx` file that includes burst onsets and offsets annotated by the human experts|
   | [lfps](https://github.com/jeelabKIST/Cho2022_BurstDetection/tree/main/data/experimental_data/ESCAPE/lfps)| Contains a `.mat` file that stores LFP signals of the robot-based escape experiment|
   | [sample_videos](https://github.com/jeelabKIST/Cho2022_BurstDetection/tree/main/data/experimental_data/ESCAPE/sample_videos)| Contains a sample video recording of the robot-based escape experiment|   

* Simulation Data

   | File           | Description                                                     |
   | :------------- | :-------------------------------------------------------------- |
   | `HM_*.mat`     | Stores the hetamps of different metrics for each algorithm      |
   | `randseed.mat` | Stores random seeds that were used to simulate synthetic signals|
   | `DC_*.mat`     | Stores the heatmaps of detection confidence scores              |   

### Utility
* `utils` include multiple sub-directories categorized by their usage. Every function script includes a description about its inputs and outputs. Please refer to each script for details.

### Tutorials

* We **_highly recommend_** reading the tutorials if the readers want to use the burst detection algorithms introduced in our paper or apply our algorithm selection rule to their own dataset.

   | File            | Description                                                                              |
   | :-------------- | :--------------------------------------------------------------------------------------- |
   | `tutorial1.mlx` | Contains guidelines for synthetic simulation, burst detection, and performance evaluation|
   | `tutorial2.mlx` | Contains guidelines for applying the algorithm selection rule                            |

### Main Scripts
* Every `.m` script in these directories starts by configuring library paths to `utils` and `data`.
    * **NOTE**: Paths in the scripts are currently set to the ones we used. To run without errors, set them to the location where your downloaded repository is at.

        ```
        util_path = genpath('PATH_TO_UTILS')
        data_path = genpath('PATH_TO_DATA')
        addpath(util_path)
        addpath(data_path)
        ```
* For `.ipynb` and `.py` scripts, you can similarly change the path by setting `input_path`, `file_path`, or `save_path` to your desired location.
* For `.r` scripts, you can configure the path by setting `data_dir` to your desired location.

### Others
* `compute_heatmap.m`

  This script is an important file that should be run prior to executing all the main scripts (except for those in `Figure1` and `Figure2`). It outputs simulation results, which are stored as a `struct` format that contains heatmaps for different metrics and algorithms. The heatmaps used in this study are already provided in the `simulation_data` directory, so you do not have to run this script unless you have a specific range of frequency band in which you want to construct the heatmaps.
  
* `save_detection_confidence.m`

  This script additionally computes heatmaps of detection confidence scores using the stored simulation results (i.e., F1-scores and temporal concurrences).

## 🎯 Requirements
The analyses and visualizations in this paper had following dependencies:

```
MATLAB 2020a (or later)
R version 4.2.2
python==3.7.4
seaborn==0.11.1
scipy==1.4.1
numpy==1.17.2
pandas==0.25.1
```

## ✏️ Citation

To cite this paper, you can use the following information:

* APA style
  
   > Cho, S., & Choi, J. H. (2023). A guide towards optimal detection of transient oscillatory bursts with unknown parameters. Journal of Neural Engineering, 20(4), 046007. https://doi.org/10.1088/1741-2552/acdffd

* Bibtex

   ```
   @article{Cho_2023,
   	doi = {10.1088/1741-2552/acdffd},
   	url = {https://doi.org/10.1088%2F1741-2552%2Facdffd},
   	year = 2023,
   	month = {jul},
   	publisher = {{IOP} Publishing},
   	volume = {20},
   	number = {4},
   	pages = {046007},
   	author = {SungJun Cho and Jee Hyun Choi},
   	title = {A guide towards optimal detection of transient oscillatory bursts with unknown parameters},
   	journal = {Journal of Neural Engineering}
   }
   ```

## 🪪 License
Copyright (c) 2022-Present [SungJun Cho](https://github.com/scho97) and [Jee Lab](https://www.jeelab.net/). `Cho2022_BurstDetection` is a free and open-source software licensed under the [MIT License](https://github.com/jeelabKIST/Cho2022_BurstDetection/blob/main/LICENSE).
