# ETH-ScopeM_Bianco
Measurement of nuclear fluo. int. of GFP-tagged proteins in *S. cerevisiae*

## Installation
1 - Download this GitHub repository  

2 - Install miniconda:  
https://docs.conda.io/en/latest/miniconda.html  

3 - Run conda prompt and install mamba in base conda environment:  
`conda install mamba -n base -c conda-forge`  

4 - Run conda prompt, from downloaded repository, and install conda environment:  
`mamba env create -f environment.yml`   

5 - Activate conda environment:  
`conda activate ETH-ScopeM_Bianco`  

6 - (Optional) Install spyder IDE:  
`pip install spyder` 

## Parameters
- `minDist` is the minimum distance between two detected nuclei (pixels)
- `minProm` is the minimum prominence (brightness) for nuclei detection (A.U.)
- `threshCoeff` nuclei segmentation threshold = `minProm` * `threshCoeff`

## Outputs
- `nData.csv`
- `nDisplay.tif`




