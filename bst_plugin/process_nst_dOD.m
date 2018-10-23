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
    sProcess.Comment     = 'raw to delta OD';
    sProcess.FileTag     = ' | dOD';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1004; %0: not shown, >0: defines place in the list of processes
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
options.option_baseline_method.Value   = {1, {'mean', 'median'}};    % {Default index, {list of entries}}

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
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
        events = sDataRaw.F.events;
    end
    
    parameters = parse_options(sProcess, sDataIn);
    
    % Remove bad channels: they won't enter dOD computation so no need to keep them     
    % Separate NIRS channels from others (NIRS_AUX etc.)
    to_keep = sDataIn.ChannelFlag ~= -1 & strcmpi({ChanneMat.Channel.Type}, 'NIRS')';
    
    % Apply dOD computation
    nirs_dOD = Compute(sDataIn.F(to_keep, :)', parameters);

    sStudy = bst_get('Study', sInputs.iStudy);
        
    % Save time-series data
    final_nirs = sDataIn.F;
    final_nirs(to_keep, :) = nirs_dOD';
    sDataOut = db_template('data');
    sDataOut.F            = final_nirs;
    sDataOut.Comment      = 'NIRS dOD';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = 1;
    sDataOut.Events       = events;
    sDataOut.DisplayUnits = 'delta OD';

    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_hb');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs.iStudy, OutputFile, sDataOut);
end


function delta_od = Compute(nirs_sig, parameters)
%% Normalize given nirs signal
% Args:
%    - nirs_sig: matrix of double, size:  nb_wavelengths x nb_samples
%        NIRS signal to normalize
%   [- parameters.method]: str , choices are: 'mean' and 'median', default is 'mean'
%        Normalization method.
%        * 'mean': divide given nirs signal by its mean, for each wavelength
%        * 'median': divide given nirs signal by its median, for each wavelength
%   [- parameters.window]: 1D array of int, of size <= nb_samples
%        Mask to be applied on the temporal axis defining the window
%        where to compute the reference signal.
%
% Output: matrix of double, size: nb_channels x nb_wavelengths
%    Normalized NIRS signal.

nb_samples = size(nirs_sig, 2);

if nargin < 2 || ~isfield(parameters, 'method')
   parameters.method = 'mean'; 
end

if nargin < 2 || ~isfield(parameters, 'baseline_window')
    parameters.baseline_window = 1:nb_samples;
end

switch parameters.method
    case 'mean'
        od_ref = mean(nirs_sig(:, parameters.baseline_window), 2);
    case 'median'
        od_ref = median(nirs_sig(:, parameters.baseline_window), 2);
end

delta_od = -log10( nirs_sig ./ repmat(od_ref, 1, nb_samples) );
end
