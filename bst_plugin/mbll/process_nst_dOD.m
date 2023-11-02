function varargout = process_nst_dOD( varargin )

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
% Authors: Thomas Vincent (2015-2016)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Raw to delta OD';
    sProcess.FileTag     = ' | dOD';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'dOD and MBLL'};
    sProcess.Index       = 1303; %0: not shown, >0: defines place in the list of processes
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'}; 
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options = get_options();
end


function options = get_options(options)

if nargin < 1
    options = struct();
end

options.option_baseline_method.Comment = 'Baseline method';
options.option_baseline_method.Type    = 'combobox';
options.option_baseline_method.Value   = {1, {'mean', 'median', 'movmean'}};    % {Default index, {list of entries}}

% === Estimation time window
options.timewindow.Comment = 'Baseline window:';
options.timewindow.Type    = 'timewindow';
options.timewindow.Value   = [];

end

function parameters = parse_options(sProcess, sDataIn)

blm_idx = sProcess.options.option_baseline_method.Value{1};
parameters.baseline_method = sProcess.options.option_baseline_method.Value{2}{blm_idx};
if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
    TimeBounds = sProcess.options.timewindow.Value{1};
else
    TimeBounds = [];
end


% Get inputs
if ~isempty(TimeBounds)
    parameters.baseline_window = panel_time('GetTimeIndices', sDataIn.Time, TimeBounds);
    if isempty(parameters.baseline_window)
        bst_report('Error', sProcess, [], 'Invalid time definition.');
        parameters = [];
        return;
    end
else
    parameters.baseline_window = [];
end

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    
    % Get baseline
    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        TimeBounds = sProcess.options.timewindow.Value{1};
    else
        TimeBounds = [];
    end
    % Comment: seconds or miliseconds
    if isempty(TimeBounds)
        Comment = [sProcess.Comment, ': All file'];
    elseif any(abs(TimeBounds) > 2)
        Comment = sprintf('%s: [%1.3fs,%1.3fs]', sProcess.Comment, TimeBounds(1), TimeBounds(2));
    else
        Comment = sprintf('%s: [%dms,%dms]', sProcess.Comment, round(TimeBounds(1)*1000), round(TimeBounds(2)*1000));
    end
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
   
    % Load channel file
    ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
    
    % Load recordings
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(1).FileName);
        events = sDataIn.Events;
        isRaw  = 0;
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
        events = sDataRaw.F.events;
        isRaw  = 1;
    end
    
    parameters = parse_options(sProcess, sDataIn);
    
    % Remove bad channels: they won't enter MBLL computation so no need to keep them 
    [good_nirs, good_channel_def] = process_nst_mbll('filter_bad_channels',sDataIn.F', ChanneMat, sDataIn.ChannelFlag);
    
    % Separate NIRS channels from others (NIRS_AUX etc.)                                                
    [fnirs, fchannel_def, nirs_other, channel_def_other] = process_nst_mbll('filter_data_by_channel_type',good_nirs, good_channel_def, 'NIRS');
    
    
    if any(any(fnirs < 0))
        msg = 'Good channels contains negative values. Consider running NISTORM -> Set bad channels';
        bst_error(msg, 'dOD quantification', 0);
    return;
    end
    
    % Apply dOD computation
    nirs_dOD = Compute(fnirs', parameters);
    
    % Re-add other channels that were not changed during MBLL
    [final_dOD, ChannelMat] = process_nst_mbll('concatenate_data',nirs_dOD', fchannel_def, nirs_other, channel_def_other);

    
    sStudy = bst_get('Study', sInputs.iStudy);
    % Create new condition because channel definition is different from original one
     cond_name = sInputs.Condition;
     if strcmp(cond_name(1:4), '@raw')
        cond_name = cond_name(5:end);
     end

     if isRaw
        iStudy = db_add_condition(sInputs.SubjectName, ['@raw', cond_name, '_dOD']);
     else
         iStudy = db_add_condition(sInputs.SubjectName, [cond_name, '_dOD']);
     end
     sStudy = bst_get('Study', iStudy);
     
%     % Save channel definition
     [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
     db_set_channel(iChannelStudy, ChannelMat, 0, 0);
        
     if ~isRaw
        % Save time-series data
        sDataOut = db_template('data');
        sDataOut.F            = final_dOD';
        sDataOut.Comment      = sDataIn.Comment;
        %sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.ChannelFlag  = ones(size(final_dOD, 2), 1);
        sDataOut.Time         = sDataIn.Time;
        sDataOut.DataType     = 'recordings'; 
        sDataOut.nAvg          = 1;
        sDataOut.Events       = events;
        sDataOut.History      = sDataIn.History;
        sDataOut              = bst_history('add', sDataOut, 'process', sProcess.Comment);
        sDataOut.DisplayUnits = 'delta OD';
        
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
         % Register in database
        db_add_data(iStudy, OutputFile, sDataOut);

     else

        % Generate a new file name in the same folder
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_0raw_OD');

        sFileIn = sDataRaw.F;
        sFileOut = out_fopen(OutputFile, 'BST-BIN',sFileIn , ChannelMat);
         % Set Output sFile structure
        sOutMat.F = sFileOut;
        sOutMat.Comment = sDataIn.Comment;
        % Save new link to raw .mat file
        bst_save(OutputFile, sOutMat, 'v6');
        % Create new channel file
        db_set_channel(iStudy, ChannelMat, 2, 0);
        % Write block
        out_fwrite(sFileOut, ChannelMat, 1, [], [], final_dOD');
        % Register in BST database
        db_add_data(iStudy, OutputFile, sOutMat);
     end

     
end


function delta_od = Compute(nirs_sig, parameters)
%% Normalize given nirs signal
% Args:
%    - nirs_sig: matrix of double, size:  nb_channels x nb_samples
%        NIRS signal to normalize
%   [- parameters.method]: str , choices are: 'mean' and 'median', default is 'mean'
%        Normalization method.
%        * 'mean': divide given nirs signal by its mean, for each wavelength
%        * 'median': divide given nirs signal by its median, for each wavelength
%   [- parameters.window]: 1D array of int, of size <= nb_samples
%        Mask to be applied on the temporal axis defining the window
%        where to compute the reference signal.
%
% Output: matrix of double, size: nb_channels x nb_samples
%    Normalized NIRS signal.

nb_samples = size(nirs_sig, 2);

if nargin < 2 || ~isfield(parameters, 'baseline_method')
   parameters.baseline_method = 'mean'; 
end

if nargin < 2 || ~isfield(parameters, 'baseline_window') || isempty(parameters.baseline_window)
    parameters.baseline_window = 1:nb_samples;
end

switch parameters.baseline_method
    case 'mean'
        od_ref = mean(nirs_sig(:, parameters.baseline_window), 2);
        od_ref = repmat(od_ref, 1, nb_samples); 
    case 'median'
        od_ref = median(nirs_sig(:, parameters.baseline_window), 2);
        od_ref = repmat(od_ref, 1, nb_samples); 
    case 'movmean'
        od_ref = movmean(nirs_sig(:, parameters.baseline_window),1000, 2);

end

delta_od = -log10( nirs_sig ./ od_ref);
end
