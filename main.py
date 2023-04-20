#%% Imports -------------------------------------------------------------------

import re
import time
import numpy as np
import pandas as pd
from skimage import io 
from pathlib import Path
from skimage.draw import disk
import matplotlib.pyplot as plt
from skimage.measure import regionprops
from skimage.feature import peak_local_max
from skimage.morphology import binary_dilation, label

#%% Initialize ----------------------------------------------------------------

# Create empty dict (imgage data = iData)
iData = {
    'name': [],
    'GFP_name': [],
    'RFP_name': [],
    'GFP_img': [],
    'RFP_img': [],
    'strain': [],
    'cond': [],
    'tp': [],
    'exp': [],
    }

# Open images and extract info
for img_path in sorted(Path('data').iterdir()): 
 
    name = img_path.stem
    
    if 'GFP' in name:
    
        # Extract name & img
        GFP_name = name
        RFP_name = name.replace('GFPdual', 'RFPtriple')
        GFP_img = io.imread(img_path) 
        RFP_img = io.imread(Path('data', RFP_name + '.tif'))
        
        # Extract info
        strain = name[0:6]    
        if 'rap' in name: cond = 'rpmc'
        if 'ut' in name: cond = 'ctrl'
        if '0min' in name: tp = 0
        if '30min' in name: tp = 30
        if '60min' in name: tp = 60
        if '120min' in name: tp = 120   
        matches = re.finditer(r'_\d_', name)
        exp = str([match.group(0) for match in matches])    
        exp = int(exp.translate(str.maketrans('', '', "[]'_")))
    
        # Append iData
        iData['name'].append(GFP_name.replace('_GFPdual', ''))
        iData['GFP_name'].append(GFP_name)
        iData['RFP_name'].append(RFP_name)
        iData['strain'].append(strain)
        iData['cond'].append(cond)
        iData['tp'].append(tp)
        iData['exp'].append(exp)
        iData['GFP_img'].append(GFP_img)
        iData['RFP_img'].append(RFP_img)

#%% 



mDist = 25
thresh = 7500

# Add keys to iData
iData['nCoords'] = []
iData['cMask'] = []
iData['nMask'] = []
iData['nLabels2D'] = []
iData['nLabels'] = []
iData['nAreas'] = []
iData['nIntRFP'] = []
iData['nIntGFP'] = []

# Segment nuclei and extract info
for i, RFP_img in enumerate(np.stack(iData['RFP_img'])):
    
    # Find nuclei positions (local maxima)
    nCoords = peak_local_max(
        RFP_img, min_distance=mDist, threshold_abs=thresh
        )
    
    # Get circle mask (centered on nuclei positions)
    cMask = np.zeros_like(RFP_img, dtype=bool)
    for coord in nCoords:
        rr, cc = disk(coord, mDist//2.5, shape=RFP_img.shape)
        cMask[rr, cc] = True
        
    # Get nuclei mask
    nMask = RFP_img > thresh/2
    nMask[cMask==False] = False
    
    # Get nuclei info
    nLabels = label(nMask)
    props = regionprops(nLabels, iData['RFP_img'][i])
    nIntRFP = np.stack([prop['intensity_mean'] for prop in props])
    props = regionprops(nLabels, iData['GFP_img'][i])
    nIntGFP = np.stack([prop['intensity_mean'] for prop in props])
    nLabels = np.stack([prop['label'] for prop in props])
    nAreas = np.stack([prop['area'] for prop in props])

    # Append iData
    iData['nCoords'].append(nCoords)
    iData['cMask'].append(cMask)
    iData['nMask'].append(nMask)
    iData['nLabels2D'].append(nLabels)
    iData['nLabels'].append(nLabels)
    iData['nAreas'].append(nAreas)
    iData['nIntRFP'].append(nIntRFP)
    iData['nIntGFP'].append(nIntGFP)
    
#%%

# Create empty DataFrame (nuclei data = nData)
headers = [
    'strain', 'cond', 'tp', 'exp',
    'xCoord', 'yCoord', 'nLabel', 'nArea', 'nIntRFP', 'nIntGFP',
    ]

nData = pd.DataFrame(columns=headers)

# Append nData
for i in range(len(iData['name'])):
    
    strain = iData['strain'][i]
    cond = iData['cond'][i]
    tp = iData['tp'][i]
    exp = iData['exp'][i]
    
    for n in range(iData['nCoords'][i].shape[0]):
        
        xCoord = iData['nCoords'][i][n][0]
        yCoord = iData['nCoords'][i][n][1]
        nLabel = iData['nLabels'][i][n]
        nArea = iData['nAreas'][i][n]
        nIntRFP = iData['nIntRFP'][i][n]
        nIntGFP = iData['nIntGFP'][i][n]
        
        nData.loc[len(nData)] = [
            strain, cond, tp, exp,
            xCoord, yCoord, nLabel, nArea, nIntRFP, nIntGFP,
            ]
        
#%%

from scipy.stats import ttest_ind, mannwhitneyu

# Create empty dict (plot data = pData)
pData = {
    'strain': [],
    'cond': [],
    'tp': [],
    'data': [],
    }

for strain in np.unique(nData['strain']):
    for cond in np.unique(nData['cond']):
        for tp in np.unique(nData['tp']):
            
            # Extract nIntRFP and nIntGFP data
            data = nData.loc[
                (nData['strain'] == strain) &
                (nData['cond'] == cond) &
                (nData['tp'] == tp),
                ['nIntRFP', 'nIntGFP']
                ]
            
            if not data.empty:
                
                pData['strain'].append(strain)
                pData['cond'].append(cond)
                pData['tp'].append(str(tp))
                pData['data'].append(data)
                
# Statistics (ttest & mwtest)
tTest_nIntRFP = []; mwTest_nIntRFP = []
tTest_nIntGFP = []; mwTest_nIntGFP = []
for strain in np.unique(nData['strain']):
    
    # Extract strain data dict (sData)
    idx = [i for i, val in enumerate(pData['strain']) if val == strain]    
    sData = {
        key: value[idx[0]:idx[-1]+1] 
        for key, value 
        in pData.items()
        }
    
    # Perform statistical tests (against ctrl 0 min)
    crtl_nIntRFP = sData['data'][0].to_numpy()[:,0]
    crtl_nIntGFP = sData['data'][0].to_numpy()[:,1]
    
    for data in sData['data']:      
        
        test_nIntRFP = data.to_numpy()[:,0]  
        _, p = ttest_ind(crtl_nIntRFP, test_nIntRFP, equal_var=True)
        tTest_nIntRFP.append(p)
        _, p = mannwhitneyu(crtl_nIntRFP, test_nIntRFP)
        mwTest_nIntRFP.append(p)
        test_nIntGFP = data.to_numpy()[:,1]
        _, p = ttest_ind(crtl_nIntGFP, test_nIntGFP, equal_var=True)
        tTest_nIntGFP.append(p)
        _, p = mannwhitneyu(crtl_nIntGFP, test_nIntGFP)
        mwTest_nIntGFP.append(p)
        
#%%
        
# Format data and labels for plotting
nIntRFP = [data.to_numpy()[:,0] for data in pData['data']]
nIntGFP = [data.to_numpy()[:,1] for data in pData['data']]

xLabels = [
    f'{tp}\n{cond}\n{strain}' 
    for (tp, cond, strain) 
    in zip(pData['tp'], pData['cond'], pData['strain'])
    ]

tLabels_nIntRFP = [
    f'p = \n{tP:.2e}' 
    for (tP, mwP) 
    in zip(tTest_nIntRFP, mwTest_nIntRFP)
    ]
tLabels_nIntRFP = [
    '' if val == 'p = \n1.00e+00' else val 
    for val in tLabels_nIntRFP 
    ]

tLabels_nIntGFP = [
    f'p = \n{tP:.2e}' 
    for (tP, mwP) 
    in zip(tTest_nIntGFP, mwTest_nIntGFP)
    ]
tLabels_nIntGFP = [
    '' if val == 'p = \n1.00e+00' else val 
    for val in tLabels_nIntGFP 
    ]
        
# Boxplot
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
box1 = ax1.boxplot(nIntRFP, labels=xLabels)
box2 = ax2.boxplot(nIntGFP, labels=xLabels)
plt.subplots_adjust(hspace=0.75)
ax1.set_title('Nuclear RFP fluo. int. (A.U.)', y=1.25)
ax2.set_title('Nuclear GFP fluo. int. (A.U.)', y=1.25)
ax2.tick_params(axis='x', labelsize=8)

# Add custom labels on top of each box
for i, (box_obj, tLabel) in enumerate(zip(box1['boxes'], tLabels_nIntRFP)):
    box_coords = box_obj.get_path().vertices
    box_x = np.mean(box_coords[:, 0])
    box_y = np.max(box_coords[:, 1])
    ax1.text(box_x, ax1.get_ylim()[1]*1.05, tLabel, ha='center', va='bottom', fontsize=8)
    
for i, (box_obj, tLabel) in enumerate(zip(box2['boxes'], tLabels_nIntGFP)):
    box_coords = box_obj.get_path().vertices
    box_x = np.mean(box_coords[:, 0])
    box_y = np.max(box_coords[:, 1])
    ax2.text(box_x, ax2.get_ylim()[1]*1.05, tLabel, ha='center', va='bottom', fontsize=8)

#%% 

# import matplotlib.pyplot as plt
# import numpy as np

# # Create some random data for the boxplot
# np.random.seed(42)
# data1 = np.random.normal(100, 10, 200)
# data2 = np.random.normal(80, 30, 200)
# data3 = np.random.normal(90, 20, 200)
# data = [data1, data2, data3]

# # Create boxplot
# fig, ax = plt.subplots()
# box = ax.boxplot(data, patch_artist=True)

# # Define custom labels
# custom_labels = ['Label 1', 'Label 2', 'Label 3']

# # Add custom labels on top of each box
# for i, (box_obj, label) in enumerate(zip(box['boxes'], custom_labels)):
#     box_coords = box_obj.get_path().vertices
#     box_x = np.mean(box_coords[:, 0])
#     box_y = np.max(box_coords[:, 1])
#     ax.text(box_x, box_y + 2, label, ha='center', va='bottom', fontsize=10)

# # Customize the plot appearance
# ax.set_xticklabels(['Data 1', 'Data 2', 'Data 3'])
# ax.set_ylabel('Values')
# ax.set_title('Boxplot with Custom Labels')

# plt.show()

#%%

# import napari
# viewer = napari.Viewer()
# viewer.add_image(np.stack(iData['GFP_img']))
# viewer.add_image(np.stack(iData['RFP_img']))
# viewer.add_image(np.stack(iData['cMask']), blending='additive', colormap='red')
# viewer.add_image(np.stack(iData['nMask']), blending='additive')
# viewer.add_image(np.stack(iData['nDisplay']), blending='additive')