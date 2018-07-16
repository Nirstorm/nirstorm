function varargout = process_nst_sci( varargin )
% process_nst_sci: compute the Scalp Coupling Index
%
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
% Authors: Thomas Vincent, 2018
%  
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute SCI';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-Scalp-coupling-index';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','raw'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    
    sProcess.options.option_low_cutoff.Comment = 'Bandpass Low cut-off';
    sProcess.options.option_low_cutoff.Type    = 'value';
    sProcess.options.option_low_cutoff.Value   = {0.5, 'Hz', 4};
    
    sProcess.options.option_high_cutoff.Comment = 'Bandpass High cut-off';
    sProcess.options.option_high_cutoff.Type    = 'value';
    sProcess.options.option_high_cutoff.Value   = {2.5, 'Hz', 4}; 
    
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInput)
    OutputFiles={};
end


