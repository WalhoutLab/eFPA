# eFPA Single-Cell Analysis

## Overview

This guide describes how to perform enhanced Flux Potential Analysis (eFPA) on specific reactions within selected individual cells using the myCobra module. The myCobra module enhances COBRApy functionalities, enabling convenient use of Flux Balance Analysis (FBA) and various data integration methods such as FPA (PMID: 33022146), eFPA (this study), Compass-(PMID: 34216539, this study), IMAT (PMID: 18711341), and IMAT++ (PMID: 33022146). It's important to note that, unlike the rest of the repository which is in MATLAB, the single-cell analysis component is implemented in Python.

## Installation

Place the Python files from this folder into your working directory and set up a Python environment with the dependencies listed in the efpa.yml file. Below is a step-by-step guide using Conda:

1. **Clone the repository:**
   
```bash 
   git clone https://github.com/YourUsername/eFPA.git
   cd eFPA/5_single-cell_analysis/
```

2. **Create a conda environment:**
   Using the efpa.yml file provided, create a Conda environment that includes all necessary dependencies.
   
```bash
   conda env create -f efpa.yml
```

3. **Activate the environment:**
   
```bash
   conda activate efpa
```

4. **Install a Solver:**
   To perform optimizations, a solver must be installed. By default, Python's GLPK library will be used if no solver is installed, but it is not optimized for our tools. We recommend using Gurobi, as it has been extensively tested with our applications. If you opt for a different advanced solver, please adjust the `<tols>` global variable in `myCobra` for optimal performance. Below is an example of how to install the Gurobi Python interface in your environment (assuming Gurobi is already installed on your computer):
   
```bash
   pip install gurobipy
```

5. **Prepare Input Data:**
   Download and extract the contents of Input.zip from [Zenodo](https://zenodo.org/records/13801228) into the Input directory. This is an example for Windows users; adjust as necessary for your OS.
   
```bash
   unzip Input.zip -d input/
```

## Usage

Execute the sc_efpa.py script from the command line to run eFPA on specified cells and reactions:

```bash
python sc_efpa.py "rxn1,rxn2,rxn3" "cell1,cell2,cell3" "output_filename"
```

- Replace "rxn1,rxn2,rxn3" with your target reactions.
- Replace "cell1,cell2,cell3" with your target cell identifiers.
- "output_filename" is where results will be saved (relative to the Output directory unless specified otherwise).

## Example run

Here is an example run and its output:

```bash
python sc_efpa.py MAR05398,MAR01380 TSP4_Muscle_diaphragm_SS2_B112866_B134045_MuscleStemCell_J19_L002,TSP4_Muscle_rectusabdominus_SS2_B114867_B134050_MuscleStemCell_J2_L004 efpa_example.pkl
```

This command analyzes specified reactions in selected muscle cells. Adjust the identifiers according to your dataset.

## Output

After running the example command, the results are saved to `"efpa_example.pkl"`. Here's how you can load and view the results:

```bash
python
import myIO
data = myIO.loadobj('Output/your_output_file.pkl')
print(data)
```

This script will load and display the results from your eFPA analysis in a Python session. Here is what the output should look like, illustrating the flux potential analysis conducted:
{'MAR05398_f': {'rFP': {'TSP4_Muscle_diaphragm_SS2_B112866_B134045_MuscleStemCell_J19_L002': 0.13324461138208257, 'TSP4_Muscle_rectusabdominus_SS2_B114867_B134050_MuscleStemCell_J2_L004': 0.4027750235549703, 'super_cond': 1.0}, 'FP': {'TSP4_Muscle_diaphragm_SS2_B112866_B134045_MuscleStemCell_J19_L002': 0.008165373471484737, 'TSP4_Muscle_rectusabdominus_SS2_B114867_B134050_MuscleStemCell_J2_L004': 0.02468248778092531, 'super_cond': 0.061281078362488556}}, 'MAR01380_f': {'rFP': {'TSP4_Muscle_diaphragm_SS2_B112866_B134045_MuscleStemCell_J19_L002': 0.37138911367066796, 'TSP4_Muscle_rectusabdominus_SS2_B114867_B134050_MuscleStemCell_J2_L004': 0.33299885201371854, 'super_cond': 1.0}, 'FP': {'TSP4_Muscle_diaphragm_SS2_B112866_B134045_MuscleStemCell_J19_L002': 0.030050699534704923, 'TSP4_Muscle_rectusabdominus_SS2_B114867_B134050_MuscleStemCell_J2_L004': 0.026944377416887808, 'super_cond': 0.08091432524151099}}}

## Documentation

For detailed usage instructions, refer to the sc_efpa.py script and the myCobra module documentation. A comprehensive tutorial on single-cell analysis will be provided soon.
