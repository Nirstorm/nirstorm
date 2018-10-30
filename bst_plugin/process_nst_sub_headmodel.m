function varargout = process_nst_sub_headmodel( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Thomas Vincent (2018)

eval(macro_method);
end


function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Cortical projection - MNE';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'NIRS - wip';
sProcess.Index       = 1206;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'results', 'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 1;


sProcess.options.head_model_fn.Comment = 'Head model file name: ';
sProcess.options.head_model_fn.Type = 'text';
sProcess.options.head_model_fn.Hidden = 1;
sProcess.options.head_model_fn.Value = '';
end