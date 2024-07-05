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