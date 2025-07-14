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
sProcess.OutputTypes = {'results'};
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
        'MriIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}

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
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

% Load subject
subjectName         = sProcess.options.subjectname.Value;
sInputs.SubjectName = subjectName;
sSubject            = bst_get('Subject', subjectName);

iStudy = db_add_condition(subjectName, sProcess.options.out_name.Value);
sStudy = bst_get('Study', iStudy);

% Load anat mri
fMRI_vol = in_mri(sProcess.options.fMRI.Value{1});
fMRI_map = fMRI_vol.Cube;

%% using Voronoi partitionning
bst_progress('start', 'Volume projection','projecting...');

SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
sCortex     = in_tess_bst(file_fullpath(sSubject.Surface(sSubject.iCortex).FileName));
nb_nodes    = size(sCortex.Vertices, 1);

voronoi         = process_nst_import_head_model('get_voronoi',sProcess, sInputs);
voronoi_mask    = (voronoi > -1) & ~isnan(voronoi);

if ~( size(fMRI_map,1) == size(voronoi,1) &&...
      size(fMRI_map,2) == size(voronoi,2) && ...
      size(fMRI_map,3) == size(voronoi,3) )

    bst_error('The dimension of the voronoi interpolator doesnt match the dimension of the projected volume');
    return;
end

func = str2func(sProcess.options.method.Value);
if ndims(fMRI_map) == 4
    surf_maps = zeros(nb_nodes,size(fMRI_map,4));
    for iTime = 1:size(fMRI_map,4)
        fmri_tmp = fMRI_map(:,:,:, iTime);

        surf_maps_tmp = accumarray(voronoi(voronoi_mask), fmri_tmp(voronoi_mask), [nb_nodes+1 1],func); 
        surf_maps_tmp(end)=[]; % trash last column

        surf_maps(:,iTime) = surf_maps_tmp;
    end

    TR = fMRI_vol.Header.dim.pixdim(5)/1000; % to check
    time = TR*(0:(size(fMRI_map,4)-1));

    sChannels = db_template('channelmat');
    sChannels.Comment = 'channel';
    sChannels.Channel = db_template('channeldesc');
    sChannels.Channel.Name = 'data';
    sChannels.Channel.Type = 'eeg';
    
    db_set_channel(iStudy, sChannels)


    sDataOut = db_template('data');
    sDataOut.F            = zeros(1, length(time));
    sDataOut.Comment      = 'BOLD data';
    sDataOut.ChannelFlag  = [1];
    sDataOut.Time         = time;
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    
    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_bold');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudy, OutputFile, sDataOut);
else
    surf_maps = accumarray(voronoi(voronoi_mask), fMRI_map(voronoi_mask), [nb_nodes+1 1],func); 
    surf_maps(end)=[]; % trash last column
    time = [1];
end

OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                                           'results_volume_projection');
[A, volName, ext] = fileparts(sProcess.options.fMRI.Value{1});

% ===== CREATE FILE STRUCTURE =====
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = volName;
ResultsMat.Function      = sProcess.options.method.Value;
ResultsMat.Time          = time;
ResultsMat.ImageGridAmp  = surf_maps;
ResultsMat.DisplayUnits  = 'BOLD';
ResultsMat.SurfaceFile   = SurfaceFile;
if ndims(fMRI_map) == 4
    ResultsMat.DataFile = sDataOut.FileName;
end

% History
ResultsMat = bst_history('add', ResultsMat, 'compute', sprintf('Projection of %s onto %s (%s)', volName,sCortex.Comment, sProcess.options.method.Value));
% Save new file structure
bst_save(OutputFile, ResultsMat, 'v6');
% Update database
db_add_data(iStudy, OutputFile, ResultsMat);
OutputFiles = {OutputFile};

bst_progress('stop');
end