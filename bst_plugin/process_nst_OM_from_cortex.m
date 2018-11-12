function varargout = process_nst_OM_from_cortex( varargin )

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
% Authors: Thomas Vincent, Alexis Machado, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute optimal montage from cortex';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 998;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

sProcess.options.cortex_to_scalp_extent.Comment = 'Extent of cortical ROI to scalp projection:';
sProcess.options.cortex_to_scalp_extent.Type = 'value';
sProcess.options.cortex_to_scalp_extent.Value = {4.0, 'cm', 2};

sProcess.options = process_nst_OM_from_head('add_OM_options', sProcess.options);
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
cortex_scout_selection = nst_get_option_selected_scout(sProcess.options, 'roi');
[head_vertex_ids, sHead] = process_nst_cpt_fluences_from_cortex('proj_cortex_scout_to_scalp', ...
                                                                 cortex_scout_selection, ...
                                                                 sProcess.options.cortex_to_scalp_extent.Value{1}.*0.01, 1);
OutputFiles = process_nst_OM_from_head('Compute', sProcess, cortex_scout_selection.sSubject, sHead, head_vertex_ids);
end
