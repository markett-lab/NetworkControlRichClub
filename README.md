# NetworkControlRichClub

This repository contains MATLAB scripts and example data used in the systematic analysis presented in the research paper "Exploring the Role of the Rich Club in Network Control of Neurocognitive States."

### Toy Example Data
Toy data from three fictive participants are included in the toy_data/ folder. These randomly generated data serve as a demonstration of the toolbox's capabilities. While the connectome sparsity is inspired by a real subject, the connectome weights are entirely random.

### Spin-Test Model-Based Null Data
The repository includes spin-test model-based null data for subset analyses. The spin_rotations.mat file contains 50 spin-rotated Lausanne ROI IDs. These rotations are based on the standard Lausanne parcellation and are suitable for general analyses.

### Example Script
The example_nct_main.m script is designed to work with the provided toy data. By adjusting the file paths in the script, you can execute all analyses with the toy data, and results will be saved automatically in the toy_data/ folder. Sample results are also included for reference.

## Prerequisites

### Required MATLAB Toolboxes and External Software

The following MATLAB toolboxes and external software are required:

- **MATLAB Toolboxes**:
  - **Signal Processing Toolbox**: Required for functions related to signal analysis and filtering.
  - **Statistics and Machine Learning Toolbox**: Used for statistical analysis, including rank-ordering and regression models.
  - **Image Processing Toolbox**: Utilized for operations involving image data or matrices that require specific image processing functions.
  - **Parallel Computing Toolbox**: Required for parallel processing to speed up computations.
  - **Optimization Toolbox**: Necessary for performing optimization tasks.

- **External Software**:
  - **GIFTI Library (for MATLAB)**: Used for reading and writing `.gii` files.
  - **Brain Connectivity Toolbox (BCT)**: A toolbox for analyzing  brain networks. The scripts use functions from BCT for graph-theoretic measures.
   - **Urs Braun's Network Control Toolbox**: Required for network control analysis. 


### Installation

1. **MATLAB Toolboxes**: Ensure that the required toolboxes are installed and licensed. You can check your installed toolboxes in MATLAB using the `ver` command.

2. **GIFTI Library**: If not already installed, you can download and install the GIFTI library for MATLAB from the [official repository](https://www.artefact.tk/software/matlab/gifti/).

3. **Brain Connectivity Toolbox (BCT)**: Download and set up BCT by following the instructions on the [official website](https://www.nitrc.org/projects/bct/).
4. 
5. **Urs Braun's Network Control Toolbox**: Clone and set up the toolbox from [this GitHub repository](https://github.com/ursbraun/network_control_and_dopamine).
