function varargout = process_nst_compute_ttest( varargin )
% process_compute_glm: compute the glm : find B such as Y = XB +e with X
% the design matrix of the experimentation. 
% This process uses OLS to fit the data.
%
% 
% 
% Further update : use more sophisticated method to fit the B
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
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
% =============================================================================@
%
% Authors: Edouard Delaire, 2018
%  
eval(macro_method);
end



%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute ttest';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'wip - GLM';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Movement_correction';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','raw'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    
    sProcess.options.fitting.Comment = 'Contrast vector';
    sProcess.options.fitting.Type    = 'combobox';
    sProcess.options.fitting.Value   = {1, {'OLS'}};
    
    

end
