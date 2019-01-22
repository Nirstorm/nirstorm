import os.path as op

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np


def plot_cmap(cm, norm, img_fn=None):
    fig = plt.figure()
    ax = fig.add_axes([0.05, 0.2, 0.9,  0.1])
    cb = mpl.colorbar.ColorbarBase(ax, cmap=cm, norm=norm,
                                   orientation='horizontal')
    cb.ax.tick_params(labelsize=38)
        
    if img_fn is not None:
        plt.savefig(img_fn, bbox_inches='tight')
        #make_white_transparent_in_img(img_fn)

def cmap_to_matlab_code(cmap, cmap_name, nb_colors=256):
    licence_comments = """% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2019 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@"""
    function_pat = 'function CMap = cmap_{}(varargin)\n\n' + \
                    licence_comments + '\n\n' + 'CMap = [\n    {}];\n'

    clist = cmap(np.linspace(0,1,nb_colors))[:,:3]
    return function_pat.format(cmap_name,
                               '\n    '.join([' '.join(['%1.5f'%c for c in rgb])
                                              for rgb in clist]))
    
cdict_hbo_np = {'red':   ((0.0, 30/255., 30/255.),
                          (0.25, 80/255., 80/255.),
                          (0.5, 125/255., 1.0),
                          (0.75, 1.0, 1.0),
                          (1.0,  1.0, 1.0)),
                'green': ((0.0,  25/255, 25/255),
                          (0.25, 0.0, 0.0),
                          (0.5,  1.0, 1.0),
                          (0.75, 190/255., 190/255.),
                          (1.0,  0.0, 0.0)),
                'blue':  ((0.0,  225/255., 225/255.),
                          (0.25, 235/255., 235/255.),
                          (0.5, 225/255., 125/255.),
                          (0.75, 0.0, 0.0),
                          (1.0,  0.0, 0.0))
}

cmap = LinearSegmentedColormap('HbOEffect', cdict_hbo_np)
vmax = 5
vmin = -5
norm = plt.Normalize(vmin, vmax) #MidpointNormalize(vmin, vmax, midpoint=0)


nb_colors = 256
cmap_name = 'ovun'
dest_dir = '/home/tom/Projects/Research/Software/Brainstorm/brainstorm3_trunk/toolbox/misc'
dest_code_fn = op.join(dest_dir, 'cmap_ovun.m')
with open(dest_code_fn, "w") as text_file:
    text_file.write(cmap_to_matlab_code(cmap, cmap_name, nb_colors))

bounds = np.linspace(0,14,15)
norm_bounds = mpl.colors.BoundaryNorm(bounds, cmap.N)
plot_cmap(cmap, norm_bounds, './cmap_ovun.png')
plt.show()
