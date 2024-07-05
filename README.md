![Python Badge](https://img.shields.io/badge/Python-3.9-rgb(69%2C132%2C182)?logo=python&logoColor=rgb(149%2C157%2C165)&labelColor=rgb(50%2C60%2C65))  
![Author Badge](https://img.shields.io/badge/Author-Benoit%20Dehapiot-blue?labelColor=rgb(50%2C60%2C65)&color=rgb(149%2C157%2C165))
![Date Badge](https://img.shields.io/badge/Created-2023--04--18-blue?labelColor=rgb(50%2C60%2C65)&color=rgb(149%2C157%2C165))
![License Badge](https://img.shields.io/badge/Licence-GNU%20General%20Public%20License%20v3.0-blue?labelColor=rgb(50%2C60%2C65)&color=rgb(149%2C157%2C165))     

# ETH-ScopeM_Bianco  
Nuclear fluo. int. GFP-tagged proteins in *S. cerevisiae*

## Index
- [Installation](#installation)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Comments](#comments)

## Installation

Pease select your operating system

<details> <summary>Windows</summary>  

### Step 1: Download this GitHub Repository 
- Click on the green `<> Code` button and download `ZIP` 
- Unzip the downloaded file to a desired location

### Step 2: Install Miniforge (Minimal Conda installer)
- Download and install [Miniforge](https://github.com/conda-forge/miniforge) for your operating system   
- Run the downloaded `.exe` file  
    - Select "Add Miniforge3 to PATH environment variable"  

### Step 3: Setup Conda 
- Open the newly installed Miniforge Prompt  
- Move to the downloaded GitHub repository
- Run the following command:  
```bash
mamba env create -f environment.yml
```
- Activate Conda environment:
```bash
conda activate Bianco
```
Your prompt should now start with `(Bianco)` instead of `(base)`

</details> 

<details> <summary>MacOS</summary>  

### Step 1: Download this GitHub Repository 
- Click on the green `<> Code` button and download `ZIP` 
- Unzip the downloaded file to a desired location

### Step 2: Install Miniforge (Minimal Conda installer)
- Download and install [Miniforge](https://github.com/conda-forge/miniforge) for your operating system   
- Open your terminal
- Move to the directory containing the Miniforge installer
- Run one of the following command:  
```bash
# Intel-Series
bash Miniforge3-MacOSX-x86_64.sh
# M-Series
bash Miniforge3-MacOSX-arm64.sh
```   

### Step 3: Setup Conda 
- Re-open your terminal 
- Move to the downloaded GitHub repository
- Run the following command: 
```bash
mamba env create -f environment.yml
```  
- Activate Conda environment:  
```bash
conda activate Bianco
```
Your prompt should now start with `(Bianco)` instead of `(base)`

</details>

## Parameters
```bash
- minDist       # minimum distance between two detected nuclei (pixels)
- minProm       # minimum prominence (brightness) for nuclei detection (A.U.)
- threshCoeff   # nuclei segmentation threshold = minProm * threshCoeff
- statTest      # statistic test should be equal to "ttest" or "utest"
```
## Outputs
```bash
- nucleiData.csv      # individual nuclei RFP and GFP fluo. int. data
- statData_tp.csv     # statistics data comparing timepoints (ttest of utest)
- statData_cond.csv   # statistics data comparing conditions (two-way ANOVA)
- statPlot.png        # plot recapitulating statistic results
- nucleiDisplay.tif   # images showing segmentation/measurments results 
    - outlines        # nuclei segmentation
    - top number      # nuclei ID
    - middle number   # RFP fluo. int.
    - bottom number   # GFP fluo. int. 
```

## Comments