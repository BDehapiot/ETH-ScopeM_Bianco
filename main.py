#%% Imports -------------------------------------------------------------------

import re
import time
from skimage import io 
from pathlib import Path

#%% Initialize ----------------------------------------------------------------

data = {
    'name': [],
    'strain': [],
    'cond': [],
    'tp': [],
    'exp': [],
    }

for img_path in sorted(Path('data').iterdir()): 
    
    # Extract name info
    name = img_path.stem
    strain = name[0:6]    
    if 'rap' in name: cond = 'rpmc'
    if 'ut' in name: cond = 'ctrl'
    if '0min' in name: tp =0
    if '30min' in name: tp = 30
    if '60min' in name: tp = 60
    if '120min' in name: tp = 120   
    matches = re.finditer(r'_\d_', name)
    exp = str([match.group(0) for match in matches])    
    exp = int(exp.translate(str.maketrans('', '', "[]'_")))
    
    # Extract img
    if 'GFP' in name:
        img_gfp = io.imread(img_path)
    if 'RFP' in name:
        img_rfp = io.imread(img_path)

    # Append data
    data['name'].append(name)
    data['strain'].append(strain)
    data['cond'].append(cond)
    data['tp'].append(time)
    data['exp'].append(exp)

#%% 
