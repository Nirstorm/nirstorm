function varargout = process_nst_cpt_fluences_from_cortex( varargin )

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
% Authors: Thomas Vincent, ZhengChen Cai (2017)
%
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute fluences from cortical scout';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1200;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

sProcess.options.cortex_to_scalp_extent.Comment = 'Extent of cortical ROI to scalp projection:';
sProcess.options.cortex_to_scalp_extent.Type = 'value';
sProcess.options.cortex_to_scalp_extent.Value = {4.0, 'cm', 1};

sProcess.options = nst_add_scout_sel_options(sProcess.options, 'roi', 'Cortical scout (target ROI):', ...
                                             'cortex', {'User scouts'}, 0);
sProcess.options = process_nst_cpt_fluences_from_head('append_mcxlab_options', sProcess.options);
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

function [head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout, extent_m, save_in_db)

if nargin < 3
    save_in_db = 0;
end

sSubject = cortex_scout.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
dis2head = pdist2(sHead.Vertices, sCortex.Vertices(cortex_scout.sScout.Vertices,:));
head_vertices = find(min(dis2head,[],2) < extent_m); 

% TODO: properly select atlas
exclude_scout = sHead.Atlas.Scouts(strcmp('FluenceExclude', {sHead.Atlas.Scouts.Label}));
if ~isempty(exclude_scout)
    head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
end

limiting_scout = sHead.Atlas.Scouts(strcmp('FluenceRegion', {sHead.Atlas.Scouts.Label}));
if ~isempty(limiting_scout)
    head_vertices = intersect(head_vertices, limiting_scout.Vertices);
end

if save_in_db && ...,
   ~any(strcmp(['From cortical ' cortex_scout.sScout.Label '(' num2str(extent_m*100) ' cm)']...,
    ,{sHead.Atlas.Scouts.Label}))
    scout_idx = size(sHead.Atlas.Scouts,2) + 1;
    sHead.Atlas.Scouts(scout_idx) = db_template('Scout');
    sHead.Atlas.Scouts(scout_idx).Vertices = head_vertices';
    sHead.Atlas.Scouts(scout_idx).Seed = head_vertices(1);
    sHead.Atlas.Scouts(scout_idx).Color = [0,0,0];
    sHead.Atlas.Scouts(scout_idx).Label = ['From cortical ' cortex_scout.sScout.Label ...
                                           '(' num2str(extent_m*100) ' cm)'];
    bst_save(file_fullpath(sSubject.Surface(sSubject.iScalp).FileName), sHead, 'v7');
    db_save();
end

end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
% Get scout vertices
cortex_scout_selection = nst_get_option_selected_scout(sProcess.options, 'roi');

[head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout_selection, ...
                                                              sProcess.options.cortex_to_scalp_extent.Value{1}.*0.01, 1);

OutputFiles = process_nst_cpt_fluences_from_head('Compute', sProcess, sSubject, sHead, head_vertices);
end
