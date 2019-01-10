function varargout = process_nst_save_matrix_csv( varargin )

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
% Authors: Thomas Vincent (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Save matrix to CSV';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1401;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    
    % File selection options
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'save', ...                           % Dialog type: {open,save}
        'Save matrix...', ...               % Window title
        'ExportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                         % Selection mode: {single,multiple}
        'files', ...                          % Selection mode: {files,dirs,files_and_dirs}
        {{'.csv'},   'ASCII: Comma-separated (*.csv)', 'ASCII-CSV'}, ... 
        'DataOut'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    sProcess.options.csv_file.Comment = 'CSV file:';
    sProcess.options.csv_file.Type    = 'filename';
    sProcess.options.csv_file.Value   = SelectOptions;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>

OutputFiles = {};

assert(length(sInput)==1);
assert(~isempty(sProcess.options.csv_file.Value{1}));

DataMat = in_bst_matrix(sInput.FileName);
[row_names, col_names] = process_nst_concat_matrices('get_axes_info', DataMat);
table_out = array2table(DataMat.Value, 'VariableNames', col_names, ...
                        'RowNames', row_names);
writetable(table_out, sProcess.options.csv_file.Value{1}, ...
           'Delimiter', ',', 'WriteRowNames', 1);
end


function nirs = Compute(channel_def, data, time, events)
if size(time, 1) == 1
    nirs.t = time';
else
    nirs.t = time;
end
if size(data, 2) ~= length(channel_def.Channel)
    data = data';
end

montage_info = nst_montage_info_from_bst_channels(channel_def.Channel);
pair_ichans = montage_info.pair_ichans;
src_coords = montage_info.src_pos;
det_coords = montage_info.det_pos;
pair_sd_indexes =  montage_info.pair_sd_indexes;

if isfield(channel_def.Nirs, 'Wavelengths')
    nirs.SD.Lambda = channel_def.Nirs.Wavelengths;
else
    nirs.SD.Lambda = channel_def.Nirs.Hb;
end
ichan = 1;
for imeasure=1:length(nirs.SD.Lambda)
    for ipair=1:size(pair_sd_indexes, 1)
        nirs.ml(ichan, 1) = pair_sd_indexes(ipair, 1);
        nirs.ml(ichan, 2) = pair_sd_indexes(ipair, 2);
        nirs.ml(ichan, 3) = 1;
        nirs.ml(ichan, 4) = imeasure;
        nirs.d(:, ichan) = data(:, pair_ichans(ipair, imeasure));
        ichan = ichan + 1;
    end
end

nirs.SD.MeasList = nirs.ml;
nirs.SD.SrcPos = src_coords;
nirs.SD.nSrcs = size(src_coords, 1);
nirs.SD.DetPos = det_coords;
nirs.SD.nDets = size(det_coords, 1);
nirs.aux = data(:, strcmp({channel_def.Channel.Type}, 'NIRS_AUX'));
nirs.events = events;
end
