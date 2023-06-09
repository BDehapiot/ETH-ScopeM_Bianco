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
- `minDist` - minimum distance between two detected nuclei (pixels)
- `minProm` - minimum prominence (brightness) for nuclei detection (A.U.)
- `threshCoeff` - nuclei segmentation threshold = `minProm` * `threshCoeff`
- `statTest` - statistic test should be equal to 'ttest' or 'utest'

## Outputs
- `nucleiData.csv` - individual nuclei RFP and GFP fluo. int. data
- `statData_tp.csv` - statistics data comparing timepoints (ttest of utest)
- `statData_cond.csv` - statistics data comparing conditions (two-way ANOVA)
- `statPlot.png` - plot recapitulating statistic results
- `nucleiDisplay.tif` - images showing segmentation/measurments results 
    - outlines = nuclei segmentation
    - top number = nuclei ID
    - middle number = RFP fluo. int.
    - bottom number = GFP fluo. int. 




