function varargout = process_nst_import_nirs_as_tstat( varargin )
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
% Authors: Thomas Vincent 2018
%  
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Import nirs as stat map';
    sProcess.Category    = 'Stat1';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'data'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;

    
    % File selection options
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'open', ...                           % Dialog type: {open,save}
        'Import NIRS...', ...               % Window title
        'ImportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                         % Selection mode: {single,multiple}
        'files', ...                          % Selection mode: {files,dirs,files_and_dirs}
        {{'.nirs'}, 'NIRS data (*.nirs)', 'NIRS-BRS'}, ... % Get all the available file formats
        'DataIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
    % Option: Event file
    sProcess.options.nirs_file.Comment = 'Nirs file: ';
    sProcess.options.nirs_file.Type    = 'filename';
    sProcess.options.nirs_file.Value   = SelectOptions;
    
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = '';

    sProcess.options.condition.Comment = 'Condition';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = 'nirs10';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment;
end

function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    sOutput = [];
        
    if isempty(sProcess.options.condition.Value)
       error('Empty condition');
    end
    
    nirs_file  = sProcess.options.nirs_file.Value{1};
    if ~exist(nirs_file, 'file')
        bst_error(['Nirs file not found: ' nirs_file]);
        return;
    end
    
    % Get subject name
    if isfield(sProcess.options, 'subjectname') && ~isempty(sProcess.options.subjectname.Value)
        SubjectName = file_standardize(sProcess.options.subjectname.Value);
    else
        SubjectName = sInputs.SubjectName;
    end
    if isempty(SubjectName)
        bst_report('Error', sProcess, [], 'Subject name is empty.');
        return;
    end
    
    % ===== GET SUBJECT =====
    % Get subject
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    if isempty(iSubject)
        bst_report('Error', sProcess, [], ['Subject "' SubjectName '" does not exist.']);
        return
    end
        
    % Output of statmap
    sOutput = db_template('statmat');
    sOutput.pmap         = p'; %Dummy value 
    sOutput.tmap         = t_stat';
    sOutput.df           = ones(nb_positions, 1) * df; %Dummy value 
    sOutput.Correction   = 'no';
    sOutput.Type = 'data';
    sOutput.ChannelFlag = ones(1,nb_positions);
    sOutput.Options.SensorTypes = 'NIRS';
    
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    sOutput.nComponents  = 1;
    
    sOutput.Comment = ['T-test : ' comment];
    
%     OutputFiles{1} = output_fn;
%     save(OutputFiles{1}, '-struct', 'sOutput');
%     db_add_data(iStudyIntra, OutputFiles{1}, sOutput);

%         % Saving the cB Matrix
% TODO: better save as either channel-space map or cortical map
%     sOutput_b = db_template('matrixmat');
%     sOutput_b.F           = con_mat';
%     sOutput_b.Comment     = [contrast ' B' ];
%     sOutput_b.Description = contrast;  
%     sOutput_b.ChannelFlag =  B.ChannelFlag;
%     sOutput_b.Time         = [1];
%     sOutput_b.DataType     = 'recordings';
%     sOutput_b.nAvg         = 1;
%     sOutput_b.DisplayUnits = 'mmol.l-1'; %TODO: check scaling
% 
%     sOutput_b = bst_history('add', sOutput_b, B.History, '');
%     sOutput_b = bst_history('add', sOutput_b, 'Subject Stat ', comment);
%     
%     OutputFiles{2} = bst_process('GetNewFilename', fileparts(sStudyIntra.FileName), 'data_cbeta_matrix');
%     save(OutputFiles{2}, '-struct', 'sOutput_b');
%     db_add_data(iStudyIntra, OutputFiles{2}, sOutput_b);
end
