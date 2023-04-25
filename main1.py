#%% Imports -------------------------------------------------------------------

import re
import cv2
import time
import numpy as np
import pandas as pd
from skimage import io 
from pathlib import Path
from skimage.draw import disk
import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from skimage.measure import regionprops
from skimage.feature import peak_local_max
from skimage.morphology import binary_dilation, label

''' There is something wrong with measured nIntRFP, not the good value '''

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
        iData['GFP_img'].append(GFP_img)
        iData['RFP_img'].append(RFP_img)
        iData['strain'].append(strain)
        iData['cond'].append(cond)
        iData['tp'].append(tp)
        iData['exp'].append(exp)

#%% 

mDist = 25
thresh = 7500

start = time.time()
print('Process')

# Add keys to iData
iData['nCoords'] = []
iData['cMask'] = []
iData['nMask'] = []
iData['nLabels2D'] = []
iData['nLabels'] = []
iData['nAreas'] = []
iData['nIntRFP'] = []
iData['nIntGFP'] = []

def process(RFP_img, GFP_img):
    
    # Find nuclei positions (local maxima)
    lmCoords = peak_local_max(
        RFP_img, min_distance=mDist, threshold_abs=thresh
        )
    
    # Get circle mask (centered on nuclei positions)
    cMask = np.zeros_like(RFP_img, dtype=bool)
    for coord in lmCoords:
        rr, cc = disk(coord, mDist//2.5, shape=RFP_img.shape)
        cMask[rr, cc] = True
        
    # Get nuclei mask & labels
    nMask = RFP_img > thresh/2
    nMask[cMask==False] = False
    
    # Get nuclei info
    nLabels2D = label(nMask)
    props = regionprops(nLabels2D, RFP_img)
    nIntRFP = np.stack([prop['intensity_mean'] for prop in props])
    props = regionprops(nLabels2D, GFP_img)
    nIntGFP = np.stack([prop['intensity_mean'] for prop in props])
    nCoords = np.stack([prop['centroid'] for prop in props]) 
    nLabels = np.stack([prop['label'] for prop in props])
    nAreas = np.stack([prop['area'] for prop in props])
    
    return nMask, nLabels2D, nIntRFP, nIntGFP, nCoords, nLabels, nAreas
    
# Run parallel
output_list = Parallel(n_jobs=-1)(
    delayed(process)(RFP_img, GFP_img)
    for (RFP_img, GFP_img) 
    in zip(
        np.stack(iData['RFP_img']), 
        np.stack(iData['GFP_img']),         
        )
    )

# nMask = [data[0] for data in output_list], axis=0).squeeze()

    # # Fill output dictionary
    # output_dict = {
    
    # # Parameters
    # 'binning': binning,
    # 'ridge_size': ridge_size,
    # 'thresh_coeff': thresh_coeff,
    
    # # Data
    # 'rsize': np.stack(
    #     [data[0] for data in output_list], axis=0).squeeze(),
    # 'ridges': np.stack(
    #     [data[1] for data in output_list], axis=0).squeeze(),
    # 'mask': np.stack(
    #     [data[2] for data in output_list], axis=0).squeeze(),

    # }

# # Segment nuclei and extract info
# for i, RFP_img in enumerate(np.stack(iData['RFP_img'])):
    
#     # Find nuclei positions (local maxima)
#     lmCoords = peak_local_max(
#         RFP_img, min_distance=mDist, threshold_abs=thresh
#         )
    
#     # Get circle mask (centered on nuclei positions)
#     cMask = np.zeros_like(RFP_img, dtype=bool)
#     for coord in lmCoords:
#         rr, cc = disk(coord, mDist//2.5, shape=RFP_img.shape)
#         cMask[rr, cc] = True
        
#     # Get nuclei mask
#     nMask = RFP_img > thresh/2
#     nMask[cMask==False] = False
    
#     # Get nuclei info
#     nLabels2D = label(nMask)
#     props = regionprops(nLabels2D, iData['RFP_img'][i])
#     nIntRFP = np.stack([prop['intensity_mean'] for prop in props])
#     props = regionprops(nLabels2D, iData['GFP_img'][i])
#     nIntGFP = np.stack([prop['intensity_mean'] for prop in props])
#     nCoords = np.stack([prop['centroid'] for prop in props]) # update nCoords (acc. to regionprops)
#     nLabels = np.stack([prop['label'] for prop in props])
#     nAreas = np.stack([prop['area'] for prop in props])

#     # Append iData
#     iData['nCoords'].append(nCoords)
#     iData['cMask'].append(cMask)
#     iData['nMask'].append(nMask)
#     iData['nLabels2D'].append(nLabels2D)
#     iData['nLabels'].append(nLabels)
#     iData['nAreas'].append(nAreas)
#     iData['nIntRFP'].append(nIntRFP)
#     iData['nIntGFP'].append(nIntGFP)
    
end = time.time()
print(f'  {(end-start):5.3f} s') 
    
#%%

# # Create empty DataFrame (nuclei data = nData)
# headers = [
#     'name', 'strain', 'cond', 'tp', 'exp',
#     'xCoord', 'yCoord', 'nLabel', 'nArea', 'nIntRFP', 'nIntGFP',
#     ]

# nData = pd.DataFrame(columns=headers)

# # Append nData
# for i in range(len(iData['name'])):
    
#     name = iData['name'][i]
#     strain = iData['strain'][i]
#     cond = iData['cond'][i]
#     tp = iData['tp'][i]
#     exp = iData['exp'][i]
    
#     for n in range(iData['nCoords'][i].shape[0]):
        
#         yCoord = iData['nCoords'][i][n][0]
#         xCoord = iData['nCoords'][i][n][1]
#         nLabel = iData['nLabels'][i][n]
#         nArea = iData['nAreas'][i][n]
#         nIntRFP = iData['nIntRFP'][i][n]
#         nIntGFP = iData['nIntGFP'][i][n]
        
#         nData.loc[len(nData)] = [
#             name, strain, cond, tp, exp,
#             yCoord, xCoord, nLabel, nArea, nIntRFP, nIntGFP,
#             ]
        
#%%

# from datetime import datetime
# date = datetime.now().strftime("%y%m%d-%H%M%S")

# # -----------------------------------------------------------------------------

# folder_name = f'{date}_analysis'
# folder_path = Path('data') / folder_name
# folder_path.mkdir(parents=True, exist_ok=False)

# # -----------------------------------------------------------------------------

# nData.to_csv(Path(folder_path) / 'nData.csv', index=False, float_format='%.3f')

# # -----------------------------------------------------------------------------

#%%

# from skimage.morphology import disk
# from skimage.exposure import rescale_intensity

# nDisplay = []
# for i in range(len(iData['name'])):
    
#     name = iData['name'][i]
#     nMask = iData['nMask'][i]
#     RFP_img = iData['RFP_img'][i]
#     GFP_img = iData['GFP_img'][i]
    
#     tmpDisplay = (binary_dilation(nMask, footprint=disk(2)) ^ nMask)
    
#     for n in range(iData['nCoords'][i].shape[0]):
        
#         yCoord = int(iData['nCoords'][i][n][0])
#         xCoord = int(iData['nCoords'][i][n][1])
#         nLabel = iData['nLabels'][i][n]
#         nIntRFP = iData['nIntRFP'][i][n]
#         nIntGFP = iData['nIntGFP'][i][n]
        
#         tmpDisplay = cv2.putText(
#             tmpDisplay.astype('uint16'), str(nLabel), (xCoord + 16, yCoord - 12),
#             cv2.FONT_HERSHEY_DUPLEX, 1, (1,1,1), 1, cv2.LINE_AA
#             ) 
        
#         tmpDisplay = cv2.putText(
#             tmpDisplay.astype('uint16'), str(int(nIntRFP)), (xCoord + 16, yCoord + 12),
#             cv2.FONT_HERSHEY_DUPLEX, 0.75, (1,1,1), 1, cv2.LINE_AA
#             )  
        
#     nDisplay.append(np.stack((RFP_img, GFP_img, tmpDisplay * 65535), axis=0))

# -----------------------------------------------------------------------------

# nDisplay = np.stack(nDisplay)
# display_name = name + '_display.tif'

# val_range = np.arange(256, dtype='uint8')
# lut_gray = np.stack([val_range, val_range, val_range])
# lut_green = np.zeros((3, 256), dtype='uint8')
# lut_green[1, :] = val_range
# lut_magenta= np.zeros((3, 256), dtype='uint8')
# lut_magenta[[0,2],:] = np.arange(256, dtype='uint8')

# io.imsave(
#     Path(folder_path, display_name),
#     nDisplay,
#     check_contrast=False,
#     imagej=True,
#     metadata={
#         'axes': 'ZCYX', 
#         'mode': 'composite',
#         'LUTs': [lut_magenta, lut_green, lut_gray],
#         }
#     )
    
#%%

# from scipy.stats import ttest_ind, mannwhitneyu

# # Create empty dict (plot data = pData)
# pData = {
#     'strain': [],
#     'cond': [],
#     'tp': [],
#     'data': [],
#     }

# for strain in np.unique(nData['strain']):
#     for cond in np.unique(nData['cond']):
#         for tp in np.unique(nData['tp']):
            
#             # Extract nIntRFP and nIntGFP data
#             data = nData.loc[
#                 (nData['strain'] == strain) &
#                 (nData['cond'] == cond) &
#                 (nData['tp'] == tp),
#                 ['nIntRFP', 'nIntGFP']
#                 ]
            
#             if not data.empty:
                
#                 pData['strain'].append(strain)
#                 pData['cond'].append(cond)
#                 pData['tp'].append(str(tp))
#                 pData['data'].append(data)
                
# # Statistics (ttest & mwtest)
# tTest_nIntRFP = []; mwTest_nIntRFP = []
# tTest_nIntGFP = []; mwTest_nIntGFP = []
# for strain in np.unique(nData['strain']):
    
#     # Extract strain data dict (sData)
#     idx = [i for i, val in enumerate(pData['strain']) if val == strain]    
#     sData = {
#         key: value[idx[0]:idx[-1]+1] 
#         for key, value 
#         in pData.items()
#         }
    
#     # Perform statistical tests (against ctrl 0 min)
#     crtl_nIntRFP = sData['data'][0].to_numpy()[:,0]
#     crtl_nIntGFP = sData['data'][0].to_numpy()[:,1]
    
#     for data in sData['data']:      
        
#         test_nIntRFP = data.to_numpy()[:,0]  
#         _, p = ttest_ind(crtl_nIntRFP, test_nIntRFP, equal_var=True)
#         tTest_nIntRFP.append(p)
#         _, p = mannwhitneyu(crtl_nIntRFP, test_nIntRFP)
#         mwTest_nIntRFP.append(p)
#         test_nIntGFP = data.to_numpy()[:,1]
#         _, p = ttest_ind(crtl_nIntGFP, test_nIntGFP, equal_var=True)
#         tTest_nIntGFP.append(p)
#         _, p = mannwhitneyu(crtl_nIntGFP, test_nIntGFP)
#         mwTest_nIntGFP.append(p)
        
#%%
        
# # Format data and labels for plotting
# nIntRFP = [data.to_numpy()[:,0] for data in pData['data']]
# nIntGFP = [data.to_numpy()[:,1] for data in pData['data']]

# xLabels = [
#     f'{tp}\n{cond}\n{strain}' 
#     for (tp, cond, strain) 
#     in zip(pData['tp'], pData['cond'], pData['strain'])
#     ]

# tLabels_nIntRFP = [
#     f'p = \n{tP:.2e}' 
#     for (tP, mwP) 
#     in zip(tTest_nIntRFP, mwTest_nIntRFP)
#     ]
# tLabels_nIntRFP = [
#     '' if val == 'p = \n1.00e+00' else val 
#     for val in tLabels_nIntRFP 
#     ]

# tLabels_nIntGFP = [
#     f'p = \n{tP:.2e}' 
#     for (tP, mwP) 
#     in zip(tTest_nIntGFP, mwTest_nIntGFP)
#     ]
# tLabels_nIntGFP = [
#     '' if val == 'p = \n1.00e+00' else val 
#     for val in tLabels_nIntGFP 
#     ]
        
# # Boxplot
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8))
# box1 = ax1.boxplot(nIntRFP, labels=xLabels)
# box2 = ax2.boxplot(nIntGFP, labels=xLabels)
# plt.subplots_adjust(hspace=0.75)
# ax1.set_title('Nuclear RFP fluo. int. (A.U.)', y=1.25)
# ax2.set_title('Nuclear GFP fluo. int. (A.U.)', y=1.25)
# ax2.tick_params(axis='x', labelsize=8)

# # Add custom labels on top of each box
# for i, (box_obj, tLabel) in enumerate(zip(box1['boxes'], tLabels_nIntRFP)):
#     box_coords = box_obj.get_path().vertices
#     box_x = np.mean(box_coords[:, 0])
#     box_y = np.max(box_coords[:, 1])
#     ax1.text(box_x, ax1.get_ylim()[1]*1.05, tLabel, ha='center', va='bottom', fontsize=8)
    
# for i, (box_obj, tLabel) in enumerate(zip(box2['boxes'], tLabels_nIntGFP)):
#     box_coords = box_obj.get_path().vertices
#     box_x = np.mean(box_coords[:, 0])
#     box_y = np.max(box_coords[:, 1])
#     ax2.text(box_x, ax2.get_ylim()[1]*1.05, tLabel, ha='center', va='bottom', fontsize=8)

#%%

# import napari
# viewer = napari.Viewer()
# viewer.add_image(np.stack(iData['GFP_img']))
# viewer.add_image(np.stack(iData['RFP_img']))
# viewer.add_image(np.stack(iData['cMask']), blending='additive', colormap='red')
# viewer.add_image(np.stack(iData['nMask']), blending='additive')
# viewer.add_image(np.stack(iData['nDisplay']), blending='additive')
