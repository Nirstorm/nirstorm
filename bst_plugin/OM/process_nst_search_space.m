function varargout = process_nst_search_space( varargin )

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
% Authors: Jean-Eudes Bornert, 2026, Edouard Delaire, 2026

eval(macro_method);
end

function sProcess = GetDescription() 
sProcess.Comment     = 'Create search space';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Optimal Montage'};
sProcess.Index       = 1100;
sProcess.Description = '';
sProcess.InputTypes  = {'import','data', 'raw'};
sProcess.OutputTypes = {'data', 'data', 'raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;
sProcess.isSeparator = 0;

sProcess.options.subjectname.Comment = 'Subject name:';
sProcess.options.subjectname.Type    = 'subjectname';
sProcess.options.subjectname.Value   = '';

sProcess.options.fluencesCond.Comment = {'panel_nst_search_space', 'Search space options:'};
sProcess.options.fluencesCond.Type    = 'editpref';
sProcess.options.fluencesCond.Value   = [];
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 

OutputFiles = {};
SubjectName = sProcess.options.subjectname.Value;
sSubject = bst_get('Subject',SubjectName);

if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
    bst_error('No available Cortex and Head surface for this subject.');
    return;
end    

sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

i_atlas = strcmp({sCortex.Atlas.Name}, sProcess.options.fluencesCond.Value.Atlas);
cortex_to_scalp_extent = sProcess.options.fluencesCond.Value.Extent;
rois = strsplit(sProcess.options.fluencesCond.Value.ROI,',');

head_vertices = [];

for i_roi = 1:length(rois)
    i_scout = strcmp({sCortex.Atlas(i_atlas).Scouts.Label}, strtrim(rois(i_roi)));

    cortex_scout            = struct();
    cortex_scout.sSubject   = sSubject;
    cortex_scout.sScout     = sCortex.Atlas(i_atlas).Scouts(i_scout);

    head_vertices = union(head_vertices,  ...
                    bst_getoutvar(1,@proj_cortex_scout_to_scalp, cortex_scout, cortex_to_scalp_extent.*0.01, 1));
end

bst_report('Info', sProcess, [], sprintf('Cortical search space projected to scalp: %d vertices created.', length(head_vertices)));

end

%% ===== UTILS =====
function [head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout, extent_m, save_in_db)

if nargin < 3
    save_in_db = 0;
end
sSubject = cortex_scout.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

dis2head = nst_pdist(sHead.Vertices, sCortex.Vertices(cortex_scout.sScout.Vertices,:));
head_vertices = find(min(dis2head,[],2) < extent_m); 

iHeadAtlas    = find(strcmp({sHead.Atlas.Name}, 'User scouts'));
exclude_scout = sHead.Atlas(iHeadAtlas).Scouts(strcmp('FluenceExclude', {sHead.Atlas(iHeadAtlas).Scouts.Label}));
if ~isempty(exclude_scout)
    head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
end

limiting_scout = sHead.Atlas(iHeadAtlas).Scouts(strcmp('FluenceRegion', {sHead.Atlas(iHeadAtlas).Scouts.Label}));
if ~isempty(limiting_scout)
    head_vertices = intersect(head_vertices, limiting_scout.Vertices);
end

if save_in_db && ...
   ~any(strcmp(['From cortical ' cortex_scout.sScout.Label '(' num2str(extent_m*100) ' cm)'], ...
    {sHead.Atlas(iHeadAtlas).Scouts.Label}))
    
    scout_idx = size(sHead.Atlas(iHeadAtlas).Scouts,2) + 1;
    sHead.Atlas(iHeadAtlas).Scouts(scout_idx) = db_template('Scout');
    sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Vertices = head_vertices';
    sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Seed = head_vertices(1);
    sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Color = [0,0,0];
    sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Label = ['From cortical ' cortex_scout.sScout.Label ...
                                           '(' num2str(extent_m*100) ' cm)'];
    
    bst_save(file_fullpath(sSubject.Surface(sSubject.iScalp).FileName), sHead, 'v7');
    db_save();
end

end
