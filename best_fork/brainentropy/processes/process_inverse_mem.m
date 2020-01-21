function varargout = process_inverse_mem( varargin )
% PROCESS_INVERSE: Compute an inverse model.

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2014 University of Southern California & McGill University
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
% Authors: Francois Tadel, 2012-2014

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'Compute sources: BEst';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 328;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Comment
    sProcess.options.comment.Comment = 'Comment: ';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = '';
    % Option: Inverse method
    sProcess.options.method.Comment = {'Maximum Entropy on the Mean (MEM)'};
    sProcess.options.method.Type    = 'radio';
    sProcess.options.method.Value   = 1;
    sProcess.options.method.Enabled = 0;
    % Options: MNE options
    sProcess.options.mem.Comment = {'panel_brainentropy', 'Source estimation options: '};
    sProcess.options.mem.Type    = 'editpref';
    sProcess.options.mem.Value   = be_main;
    % Option: Sensors selection
    sProcess.options.sensortypes.Comment = 'Sensor types:&nbsp;&nbsp;&nbsp;&nbsp;';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'MEG, MEG MAG, MEG GRAD, EEG';
    % Option: Output
    sProcess.options.sep3.Type      = 'separator';
    sProcess.options.output.Comment = {'Full results: one per file'};
    sProcess.options.output.Type    = 'radio';
    sProcess.options.output.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    
    % ===== GET OPTIONS =====
    % Default inverse options
    OPTIONS = process_inverse('Compute');
    % MEM options
    if ~isfield(sProcess.options.mem.Value, 'MEMpaneloptions')
        fprintf('\n\n***\tError in BEst process\t***\n\tyou MUST edit options before lauching the MEM.\n\n')
        return
    end
    OPTIONS.MEMpaneloptions =   sProcess.options.mem.Value.MEMpaneloptions;
    % Get options
    OPTIONS.InverseMethod   =   'mem';
	OPTIONS.SourceOrient    =   {'fixed'};
    % Output
	iStudies = [sInputs.iStudy];
    iDatas   = [sInputs.iItem];
	OPTIONS.ComputeKernel = 0;
    % Get modalities in channel files
    AllSensorTypes = unique(cat(2, sInputs.ChannelTypes));
    AllSensorTypes = intersect(AllSensorTypes, {'MEG MAG', 'MEG GRAD', 'MEG', 'EEG', 'ECOG', 'SEEG'});
    if any(ismember(AllSensorTypes, {'MEG MAG', 'MEG GRAD'}))
        AllSensorTypes = setdiff(AllSensorTypes, 'MEG');
    end
    % Get valid modalities in head models
    allChanFiles = unique({sInputs.ChannelFile});
    for i = 1:length(allChanFiles)
        % Get study
        sStudy = bst_get('ChannelFile', allChanFiles{i});
        % Check if all the files exist
        if isempty(sStudy.Channel) || isempty(sStudy.HeadModel) || isempty(sStudy.iHeadModel) || isempty(sStudy.NoiseCov)
            bst_report('Error', sProcess, [], 'No channel file, noise covariance, or headmodel or for at least one of the files.');
            return;
        end
        % Remove all the modalities that do not exist in the headmodels
        if isempty(sStudy.HeadModel(sStudy.iHeadModel).MEGMethod)
            AllSensorTypes = setdiff(AllSensorTypes, {'MEG', 'MEG MAG', 'MEG GRAD'});
        end
        if isempty(sStudy.HeadModel(sStudy.iHeadModel).EEGMethod)
            AllSensorTypes = setdiff(AllSensorTypes, {'EEG'});
         end
        if isempty(sStudy.HeadModel(sStudy.iHeadModel).ECOGMethod)
            AllSensorTypes = setdiff(AllSensorTypes, {'ECOG'});
        end
        if isempty(sStudy.HeadModel(sStudy.iHeadModel).SEEGMethod)
            AllSensorTypes = setdiff(AllSensorTypes, {'SEEG'});
        end
    end
    % Selected sensor types
    OPTIONS.DataTypes = strtrim(str_split(sProcess.options.sensortypes.Value, ',;'));
    if ismember('MEG', OPTIONS.DataTypes) && any(ismember({'MEG GRAD','MEG MAG'}, AllSensorTypes))
        OPTIONS.DataTypes = union(setdiff(OPTIONS.DataTypes, 'MEG'), {'MEG MAG', 'MEG GRAD'});
    end
    OPTIONS.DataTypes = intersect(OPTIONS.DataTypes, AllSensorTypes);
    if isempty(OPTIONS.DataTypes)
        strTypes = '';
        for i = 1:length(AllSensorTypes)
            if (i > 1)
                strTypes = [strTypes, ', '];
            end
            strTypes = [strTypes, AllSensorTypes{i}];
        end
        bst_report('Error', sProcess, [], ['No valid sensor type selected.' 10 'Valid options are: ' strTypes]);
        return;
    end
    % Comment
    if isfield(sProcess.options, 'comment') && isfield(sProcess.options.comment, 'Value') && ~isempty(sProcess.options.comment.Value)
        OPTIONS.Comment = sProcess.options.comment.Value;
    end
    % No messages
    OPTIONS.DisplayMessages = 0;

    % ===== START COMPUTATION =====
    % Call head modeler
    [AllFiles, errMessage] = process_inverse('Compute', iStudies, iDatas, OPTIONS);
    % Report errors
    if isempty(AllFiles) && ~isempty(errMessage)
        bst_report('Error', sProcess, sInputs, errMessage);
        return;
    elseif ~isempty(errMessage)
        bst_report('Warning', sProcess, sInputs, errMessage);
    end
    % For shared kernels: Return only the source files corresponding to the recordings that were in input
    if isempty(iDatas)
        % Loop on the output files (all links): find the ones that match the input
        for iFile = 1:length(AllFiles)
            % Resolve link: get data file
            [ResFile, DataFile] = file_resolve_link(AllFiles{iFile});
            % Find one that matches the inputs
            iInput = find(file_compare({sInputs.FileName}, file_short(DataFile)));
            % If founf: add to the output files
            if ~isempty(iInput)
                OutputFiles{end+1} = AllFiles{iFile};
            end
        end
    else
        OutputFiles = AllFiles;
    end
end

            

