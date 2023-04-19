#%% Imports -------------------------------------------------------------------

import re
import time
import numpy as np
import pandas as pd
from skimage import io 
from pathlib import Path

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

from skimage.draw import disk
from skimage.feature import peak_local_max
from skimage.morphology import binary_dilation, label
from skimage.measure import regionprops

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
iData['nDisplay'] = []

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

    # Make displays
    nDisplay = binary_dilation(nMask) ^ nMask
    
    # Append iData
    iData['nCoords'].append(nCoords)
    iData['cMask'].append(cMask)
    iData['nMask'].append(nMask)
    iData['nLabels2D'].append(nLabels)
    iData['nLabels'].append(nLabels)
    iData['nAreas'].append(nAreas)
    iData['nIntRFP'].append(nIntRFP)
    iData['nIntGFP'].append(nIntGFP)
    iData['nDisplay'].append(nDisplay)
    
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

import matplotlib.pyplot as plt

# Create empty dict (tag data = tData)
tData = {
    'tag': [],
    'data': [],
    }

for strain in np.unique(nData['strain']):
    for cond in np.unique(nData['cond']):
        for tp in np.unique(nData['tp']):
            
            data = nData.loc[
                (nData['strain'] == strain) &
                (nData['cond'] == cond) &
                (nData['tp'] == tp),
                ['nIntRFP', 'nIntGFP']
                ]
            
            if not data.empty:
                
                tData['tag'].append(f'{strain}_{cond}_{tp}')
                tData['data'].append(data)
                
# -----------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axisartist.parasite_axes import SubplotHost

tags = [tag for tag in tData['tag']]
nIntRFP = [data.to_numpy()[:,0] for data in tData['data']]
nIntGFP = [data.to_numpy()[:,1] for data in tData['data']]

fig = plt.figure()
ax1 = SubplotHost(fig, 111)
fig.add_subplot(ax1)

ax1.boxplot(nIntRFP, labels=['0', '30', '60', '120', '0', '30', '60', '120'])

ax2 = ax1.twiny()
offset = 0, -25 # Position of the second axis
new_axisline = ax2.get_grid_helper().new_fixed_axis
ax2.axis["bottom"] = new_axisline(loc="bottom", axes=ax2, offset=offset)
# ax2.axis["top"].set_visible(False)
            
# tags = [tag for tag in tData['tag']]
# nIntRFP = [data.to_numpy()[:,0] for data in tData['data']]
# nIntGFP = [data.to_numpy()[:,1] for data in tData['data']]

# fig, ax = plt.subplots(nrows=2, ncols=1)
# ax[0].boxplot(nIntRFP, labels=tags)
# ax[1].boxplot(nIntGFP, labels=tags)

# for ax in fig.get_axes():
#     ax.label_outer()
    
# ax.set_xticklabels([])

#%%

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
# from mpl_toolkits.axisartist.parasite_axes import SubplotHost

# fig1 = plt.figure()
# ax1 = SubplotHost(fig1, 111)
# fig1.add_subplot(ax1)

# # Some data
# x = np.arange(1,6)
# y = np.random.random(len(x))

# # First X-axis
# ax1.plot(x, y)
# ax1.set_xticks(x)
# ax1.set_xticklabels(['dog', 'cat', 'horse', 'lizard', 'crocodile'])
# #ax1.xaxis.set_label_text('First X-axis') # Uncomment to label axis
# ax1.yaxis.set_label_text("Sample data")

# # Second X-axis
# ax2 = ax1.twiny()
# offset = 0, -25 # Position of the second axis
# new_axisline = ax2.get_grid_helper().new_fixed_axis
# ax2.axis["bottom"] = new_axisline(loc="bottom", axes=ax2, offset=offset)
# ax2.axis["top"].set_visible(False)

# ax2.set_xticks([0.0, 0.6, 1.0])
# ax2.xaxis.set_major_formatter(ticker.NullFormatter())
# ax2.xaxis.set_minor_locator(ticker.FixedLocator([0.3, 0.8]))
# ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['mammal', 'reptiles']))

# # Third X-axis
# ax3 = ax1.twiny()
# offset = 0, -50
# new_axisline = ax3.get_grid_helper().new_fixed_axis
# ax3.axis["bottom"] = new_axisline(loc="bottom", axes=ax3, offset=offset)
# ax3.axis["top"].set_visible(False)

# ax3.set_xticks([0.0, 1.0])
# ax3.xaxis.set_major_formatter(ticker.NullFormatter())
# ax3.xaxis.set_minor_locator(ticker.FixedLocator([0.5]))
# ax3.xaxis.set_minor_formatter(ticker.FixedFormatter(['vertebrates']))

# ax1.grid(1)
# plt.show()



#%%

# import napari
# viewer = napari.Viewer()
# viewer.add_image(np.stack(iData['GFP_img']))
# viewer.add_image(np.stack(iData['RFP_img']))
# viewer.add_image(np.stack(iData['cMask']), blending='additive', colormap='red')
# viewer.add_image(np.stack(iData['nMask']), blending='additive')
# viewer.add_image(np.stack(iData['nDisplay']), blending='additive')