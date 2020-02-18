function varargout = process_nst_OM_from_cap( varargin )

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
% Authors: Alexis Machado, Thomas Vincent, ZhengChen Cai (2017-2018)

% TODO: expose projection distance from cortex ROI to head
% TODO: only one wavelength (quickly test if montage change from one wl to
%      another)
% TODO: test limitation of model size -> check with 32x32 or 16x16.
% TODO: document computation time
% TODO: doc params:
%   - nb sources/det is the nb of optodes placed by algo then post proc to
%     remove positions
%   - nb adjacent is MIN (can be higher). Has to be associated to dist MAX
% TODO: fix paper reference + ref to equations

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute optimal montage from cap';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1409;
sProcess.Description = '';
sProcess.isSeparator = 1;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;




%sProcess.options = nst_add_scout_sel_options(struct(), 'head', 'Head scout (search space):', ...
%    'scalp', {'User scouts'}, 0);
sProcess.options.subject.Comment = 'Select Subject:';
sProcess.options.subject.Type    = 'subjectname';
sProcess.options.subject.Value   = [];

  
sProcess.options = add_OM_options(sProcess.options);
end

function options = add_OM_options(options)

% Add selector of cortical scout for target VOI
options = nst_add_scout_sel_options(options, 'roi', 'Cortical scout (target ROI):', ...
    'cortex', {'User scouts'}, 1);

options.condition_name.Comment = 'Output condition name:';
options.condition_name.Type = 'text';
options.condition_name.Value = '';

options.wavelengths.Comment = 'Wavelengths (nm) [coma-separated list]';
options.wavelengths.Type    = 'text';
options.wavelengths.Value = '';

options.data_source.Comment = 'Fluence Data Source (URL or path)';
options.data_source.Type    = 'text';
options.data_source.Value = [nst_get_repository_url() '/fluence/'];

options.nb_sources.Comment = 'Number of sources:';
options.nb_sources.Type = 'value';
options.nb_sources.Value = {4,'',0};

options.nb_detectors.Comment = 'Number of detectors:';
options.nb_detectors.Type = 'value';
options.nb_detectors.Value =  {8,'',0};

options.nAdjacentDet.Comment = 'Number of Adjacent:';
options.nAdjacentDet.Type = 'value';
options.nAdjacentDet.Value = {2,'',0};

options.sep_optode.Comment = 'Range of optodes distance:';
options.sep_optode.Type = 'range';
options.sep_optode.Value = {[15 55],'mm',0};



% WIP: add option to cache weight table and avoid recomputation TODO: make
% more robust and clearer usage
options.exist_weight.Hidden = 1;
options.exist_weight.Comment = 'Use existing weight tables (speed up)';
options.exist_weight.Type = 'checkbox';
options.exist_weight.Value = 0;

SelectOptions = {...
    '', ...                            % Filename
    '', ...                            % FileFormat
    'save', ...                        % Dialog type: {open,save}
    'Select output folder...', ...     % Window title
    'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
    'single', ...                      % Selection mode: {single,multiple}
    'dirs', ...                        % Selection mode: {files,dirs,files_and_dirs}
    {{'.folder'}, '*.*'}, ... % Available file formats
    'MriOut'};                         % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
% Option definition
% TODO: add flag to enable ouput
options.outputdir.Hidden = 1;
options.outputdir.Comment = 'Folder for weight table:';
options.outputdir.Type    = 'filename';
options.outputdir.Value   = SelectOptions;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

sInputs.SubjectName=sProcess.options.subject.Value;
% Get available position for optodes 
sStudy=bst_get('Study');
ChannelMat = in_bst_channel(sStudy.Channel.FileName);
sProcess.additional_channel=ChannelMat;

n_chan=size(ChannelMat.Channel,2);
loc=zeros(n_chan,3);

for i=1:n_chan
   loc(i,:)= ChannelMat.Channel(i).Loc;
   name{i}=ChannelMat.Channel(i).Name;
end    

sProcess.reference.loc=loc;
sProcess.reference.name=name;

% Load head mesh

[sSubject, iSubject] = bst_get('Subject', sInputs.SubjectName);
head_mesh_fn = sSubject.Surface(sSubject.iScalp).FileName;
sHead = in_tess_bst(head_mesh_fn);

% Load anat mri
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

% Find closest head vertices (for which we have fluence data)
% Put everything in mri referential
head_vertices_mri = cs_convert(sMri, 'scs', 'mri', sHead.Vertices) * 1000;
src_locs_mri = cs_convert(sMri, 'scs', 'mri', loc) * 1000;

%head_normals = process_nst_ComputeFluencesforOptodes('computeVertNorm',head_vertices_mri,sHead.Faces);
head_vertices = knnsearch(head_vertices_mri, src_locs_mri);



% Get scout vertices & load head mesh
%head_scout_selection = nst_get_option_selected_scout(sProcess.options, 'head');
%head_vertices = head_scout_selection.sScout.Vertices;

%sSubject = head_scout_selection.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);

OutputFiles = process_nst_OM_from_head('Compute',sProcess,sSubject, sHead, head_vertices);

end
