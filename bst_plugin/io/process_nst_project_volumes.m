function varargout = process_nst_project_volumes( varargin )

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
% Authors:  Edouard Delaire (2022)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Project Volume to Surface';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1401;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;
sProcess.isSeparator = 1;


sProcess.options.subjectname.Comment = 'Subject name';
sProcess.options.subjectname.Type    = 'subjectname';
sProcess.options.subjectname.Value   = '';

sProcess.options.out_name.Comment = 'Output Folder name';
sProcess.options.out_name.Type='text';
sProcess.options.out_name.Value = '';

SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'open', ...                        % Dialog type: {open,save}
        'Import MRI...', ...               % Window title
        'ImportAnat', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'files', ...                       % Selection mode: {files,dirs,files_and_dirs}
        bst_get('FileFilters', 'mri'), ... % Get all the available file formats
        'MriIn'};                        % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    % Option definition
    % TODO: add flag to enable ouput
sProcess.options.fMRI.Comment = 'Select Volume:';
sProcess.options.fMRI.Type    = 'filename';
sProcess.options.fMRI.Value   = SelectOptions;
    

sProcess.options.label.Comment = ['<b>Projection Function; </b> <BR>' ...
                                   'Function to apply within each Voronoi cell'];
sProcess.options.label.Type = 'label';

sProcess.options.method.Comment = {'mean','median', 'mode', 'min','max'; ...
                                   'mean','median', 'mode','min' 'max'};
sProcess.options.method.Type    = 'radio_label';
sProcess.options.method.Value   = 'mean';
    
% Smoothing options
sProcess.options.label2.Comment = ['<b>Smoothing: </b> <BR>'];
sProcess.options.label2.Type = 'label';

sProcess.options.do_smoothing.Comment = 'Apply smoothing after the projection';
sProcess.options.do_smoothing.Type    = 'checkbox';
sProcess.options.do_smoothing.Value   = 1;
sProcess.options.do_smoothing.Controller = 'smoothing';

% === DESCRIPTION   copied from process_ssmooth_surfstat
sProcess.options.help.Comment = ['This process uses SurfStatSmooth (SurfStat, KJ Worsley).<BR><BR>' ...
    'The smoothing is based only on the surface topography, <BR>' ...
    'not the real geodesic distance between two vertices.<BR>', ...
    'The input in mm is converted to a number of edges based<BR>', ...
    'on the average distance between two vertices in the surface.<BR><BR>'];
sProcess.options.help.Type    = 'label';
sProcess.options.help.Class = 'smoothing';

% === FWHM (kernel size)
sProcess.options.fwhm.Comment = '<B>FWHM</B> (Full width at half maximum):  ';
sProcess.options.fwhm.Type    = 'value';
sProcess.options.fwhm.Value   = {5, 'mm', 2};
sProcess.options.do_smoothing.Class = 'smoothing';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

do_smoothing    = sProcess.options.do_smoothing.Value;
subjectName     = sProcess.options.subjectname.Value;

% Load head mesh
[sSubject, iSubject] = bst_get('Subject', subjectName);

iStudy = db_add_condition(subjectName, sProcess.options.out_name.Value);
sStudy = bst_get('Study', iStudy);

% Load anat mri
fMRI_vol = in_mri(sProcess.options.fMRI.Value{1});
fMRI_map = fMRI_vol.Cube;

%% using Voronoi partitionning
bst_progress('start', 'Volume projection','projecting...');

SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
sCortex     = load(file_fullpath(sSubject.Surface(sSubject.iCortex).FileName));
nb_nodes    = size(sCortex.Vertices, 1);

voronoi         = get_voronoi(subjectName);
voronoi_mask    = (voronoi > -1) & ~isnan(voronoi);

if ~( size(fMRI_map,1) == size(voronoi,1) &&...
      size(fMRI_map,2) == size(voronoi,2) && ...
      size(fMRI_map,3) == size(voronoi,3) )

    bst_error('The dimension of the voronoi interpolator doesnt match the dimension of the projected volume');
    return;
end

func = str2func(sProcess.options.method.Value);
if ndims(fMRI_map) == 4
    surf_maps = zeros(size(nb_nodes),size(fMRI_map,4));
    for iTime = 1:length(size(fMRI_map,4))
        fmri_tmp = fMRI_map(:,:,:, iTime);

        surf_maps_tmp = accumarray(voronoi(voronoi_mask), fmri_tmp(voronoi_mask), [nb_nodes+1 1],func); 
        surf_maps_tmp(end)=[]; % trash last column

        surf_maps(:,iTime) = surf_maps_tmp;
    end

    TR = fMRI_vol.Header.dim.pixdim(4); % to check
    time = TR*(0:size(fMRI_map,4));
else
    surf_maps = accumarray(voronoi(voronoi_mask), fMRI_map(voronoi_mask), [nb_nodes+1 1],func); 
    surf_maps(end)=[]; % trash last column
    time = [1];
end

[sSubject, iSubject] = bst_get('Subject', subjectName);

if do_smoothing
    FWHM = sProcess.options.fwhm.Value{1} / 1000;
    cortex_mesh = sSubject.Surface(sSubject.iCortex).FileName;
    sCortex = in_tess_bst(cortex_mesh);
    surf_maps = surface_smooth(FWHM, sCortex, surf_maps,1);
else
    FWHM = 0;
end

[A,B] = fileparts(sProcess.options.fMRI.Value{1});
[sStudy, ResultFile] = add_surf_data_fMRI(surf_maps, time, ...
    SurfaceFile, sprintf('Projection %s (%s)', B, sProcess.options.method.Value), ...
    iStudy, sStudy,  ...
    sprintf('Projection %s (%s - FWHM %d)', B, sProcess.options.method.Value,FWHM));

OutputFiles{end+1} = ResultFile;
bst_progress('stop');

end

function [sStudy, ResultFile] = add_surf_data_fMRI(data, time, SurfaceFile, name, ...
    iStudy, sStudy, history_comment)

%% Save a cortical map to brainstorm with given data

ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
    ['results_' protect_fn_str(name)]);

% ===== CREATE FILE STRUCTURE =====
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = name;
ResultsMat.Function      = '';
ResultsMat.ImageGridAmp = data;
ResultsMat.Time          = time;
ResultsMat.DataFile      = [];
% ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
% ResultsMat.HeadModelType = head_model.HeadModelType;
ResultsMat.ChannelFlag   = [];
ResultsMat.GoodChannel   = [];
ResultsMat.SurfaceFile   = file_short(SurfaceFile);
ResultsMat.GridLoc    = [];
ResultsMat.GridOrient = [];
ResultsMat.nAvg      = 1;
% History
ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
% Save new file structure
bst_save(ResultFile, ResultsMat, 'v6');
% ===== REGISTER NEW FILE =====
% Create new results structure
newResult = db_template('results');
newResult.Comment       = name;
newResult.FileName      = file_short(ResultFile);
newResult.DataFile      = ''; %sInputs.FileName;
newResult.isLink        = 0;
% newResult.HeadModelType = ResultsMat.HeadModelType;
% Add new entry to the database
iResult = length(sStudy.Result) + 1;
sStudy.Result(iResult) = newResult;
% Update Brainstorm database
bst_set('Study', iStudy, sStudy);
end

function voronoi = get_voronoi(SubjectName)
[sSubject, iSubject] = bst_get('Subject', SubjectName);
voronoi_fn = process_mfip_compute_voronoi('get_voronoi_fn', sSubject);

voronoi_bst = in_mri_bst(voronoi_fn);
voronoi = voronoi_bst.Cube;
end

function fn = protect_fn_str(sfn)
fn = strrep(sfn, ' ', '_');
fn = strrep(fn, '"', '');
fn = strrep(fn, ':', '_');
fn = strrep(fn, '(', '_');
fn = strrep(fn, ')', '_');
end

%% ===== surface smoothing  =====
function sens_smoothed = surface_smooth(FWHM, SurfaceMat, sens_temp, dispInfo)
% ===== PROCESS =====
% Convert surface to SurfStat format
cortS.tri = SurfaceMat.Faces;
cortS.coord = SurfaceMat.Vertices';

% Get the average edge length
[vi,vj] = find(SurfaceMat.VertConn);
Vertices = SurfaceMat.VertConn;
meanDist = mean(sqrt((Vertices(vi,1) - Vertices(vj,1)).^2 + (Vertices(vi,2) - Vertices(vj,2)).^2 + (Vertices(vi,3) - Vertices(vj,3)).^2));
% FWHM in surfstat is in mesh units: Convert from millimeters to "edges"
FWHMedge = FWHM ./ meanDist;

% Display the result of this conversion
msgInfo = ['Average distance between two vertices: ' num2str(round(meanDist*10000)/10) ' mm' 10 ...
    'SurfStatSmooth called with FWHM=' num2str(round(FWHMedge * 1000)/1000) ' edges'];

%bst_report('Info', sProcess, sInput, msgInfo);
if dispInfo
    disp(msgInfo);
    disp(['SMOOTH> ' strrep(msgInfo, char(10), [10 'SMOOTH> '])]);
end

if round(FWHMedge * 1000)/1000 ==0
    disp('WARNING: FWHM too small, smoothing will not be performed.')
end
% Smooth surface

sens_smoothed = SurfStatSmooth(sens_temp', cortS, FWHMedge)';

end
