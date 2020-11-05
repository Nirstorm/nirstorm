function varargout = process_nst_glm_fit1( varargin )
% process_compute_glm: compute the glm : find B such as Y = XB +e with X
% 
% OlS_fit use an ordinary least square algorithm to find B : B= ( X^{T}X)^{-1} X^{T} Y 
% AR-IRLS : Details about AR-IRLS algorithm can be found here :
%   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3756568/ 
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
% Authors: Edouard Delaire, Thomas Vincent 2018-2020
%
% TODO: use different output types for channel-space and surface analyses
%       -> more convenient for contrast computation afterwards, to be able to map
%       input to output types
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>


    sProcess = process_nst_glm_fit('GetDescription');
    
    sProcess.Category    = 'File';
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = process_nst_glm_fit('Run',sProcess,sInput,[]);
end
