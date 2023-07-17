function varargout = process_nst_mbll_dOD( varargin )

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
% Authors: Edouard Delaire(2020)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'MBLL - delta OD to delta [HbO], [HbR] & [HbT]';
    sProcess.FileTag     = '_Hb';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'dOD and MBLL'};
    sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Compute_.5BHb.5D_variations_-_Modified_Beer-Lambert_Law';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'}; %TODO: 'raw' -> 'raw' or 'raw' -> 'data'?
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    sProcess.options.option_age.Comment = 'Age';
    sProcess.options.option_age.Type    = 'value';
    sProcess.options.option_age.Value   = {25, 'years', 2};
    
    sProcess.options.option_pvf.Comment = 'PVF';
    sProcess.options.option_pvf.Type    = 'value';
    sProcess.options.option_pvf.Value   = {50, '', 0};
    
    sProcess.options.option_do_plp_corr.Comment = 'DPF correction';
    sProcess.options.option_do_plp_corr.Type    = 'checkbox';
    sProcess.options.option_do_plp_corr.Value   = 1;
    
    dpf_method_choices = process_nst_mbll('get_dpf_method_choices');
    
    sProcess.options.option_dpf_method.Comment = 'DPF method';
    sProcess.options.option_dpf_method.Type    = 'combobox';
    sProcess.options.option_dpf_method.Value   = {dpf_method_choices.DUNCAN1996, ...
                                                  fieldnames(dpf_method_choices)};    % {Default index, {list of entries}}
    
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment =  process_nst_dOD('FormatComment', sProcess);
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFile = {};
    
    % TODO: check for negative values
    % Get option values   
    age             = sProcess.options.option_age.Value{1};

    do_plp_corr     = sProcess.options.option_do_plp_corr.Value;
    pvf = sProcess.options.option_pvf.Value{1};
    dpf_method = sProcess.options.option_dpf_method.Value{1};

    % Load channel file
    ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
    
    % Load recordings
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(1).FileName);
        events = sDataIn.Events;
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
        events = sDataRaw.F.events;
    end
    

    
    % Remove bad channels: they won't enter MBLL computation so no need to keep them 
    [good_nirs, good_channel_def] = process_nst_mbll('filter_bad_channels',sDataIn.F', ChanneMat, sDataIn.ChannelFlag);
    
    % Separate NIRS channels from others (NIRS_AUX etc.)                                                
    [fnirs, fchannel_def, nirs_other, channel_def_other] = process_nst_mbll('filter_data_by_channel_type',good_nirs, good_channel_def, 'NIRS');
    

    % Apply MBLL
    % TODO: add baseline window and expose it
    [nirs_hb, channels_hb] = process_nst_mbll('Compute',fnirs, fchannel_def, age, [], do_plp_corr, pvf, dpf_method); 
    
    % Re-add other channels that were not changed during MBLL
    [final_nirs, ChannelMat] = process_nst_mbll('concatenate_data',nirs_hb, channels_hb, nirs_other, channel_def_other);

    if ~isempty(sDataIn.Std)
         % Remove bad channels: they won't enter MBLL computation so no need to keep them 
        [good_nirs, good_channel_def] = process_nst_mbll('filter_bad_channels',sDataIn.Std', ChanneMat, sDataIn.ChannelFlag);
        
        % Separate NIRS channels from others (NIRS_AUX etc.)                                                
        [fnirs, fchannel_def, nirs_other, channel_def_other] = process_nst_mbll('filter_data_by_channel_type',good_nirs, good_channel_def, 'NIRS');
        
    
        % Apply MBLL
        % TODO: add baseline window and expose it
        [nirs_hb, channels_hb] = process_nst_mbll('Compute',fnirs, fchannel_def, age, [], do_plp_corr, pvf, dpf_method); 
        
        % Re-add other channels that were not changed during MBLL
        [final_nirs_std, ~] = process_nst_mbll('concatenate_data',nirs_hb, channels_hb, nirs_other, channel_def_other);
    
    end

    % Create new condition because channel definition is different from original one
    cond_name = sInputs.Condition;
    if length(cond_name)>=4 && strcmp(cond_name(1:4), '@raw')
        cond_name = cond_name(5:end);
    end
    iStudy = db_add_condition(sInputs.SubjectName, [cond_name, '_Hb']);
    sStudy = bst_get('Study', iStudy);
    
    % Save channel definition
    [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
    db_set_channel(iChannelStudy, ChannelMat, 0, 0);
    
    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = final_nirs'; % TOCHECK brainstorm expects MILLIMOL!!
    if ~isempty(sDataIn.Std)
        sDataOut.Std          = final_nirs_std';
    end
    sDataOut.Comment      = sDataIn.Comment;
    sDataOut.ChannelFlag  = ones(size(final_nirs, 2), 1);
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = sDataIn.nAvg;
    sDataOut.Leff         = sDataIn.Leff;
    sDataOut.Events       = events;
    sDataOut.History      = sDataIn.History;
    sDataOut = bst_history('add', sDataOut, 'process', sProcess.Comment);
    sDataOut.DisplayUnits = 'mol.l-1';

    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_hb');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudy, OutputFile, sDataOut);
end

