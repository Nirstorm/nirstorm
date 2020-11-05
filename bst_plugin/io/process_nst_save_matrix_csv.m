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
    sProcess.Index       = 1901;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    
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
    
    sProcess.options.ignore_cols_all_zeros.Comment = 'Ignore all-zero columns';
    sProcess.options.ignore_cols_all_zeros.Type    = 'checkbox';
    sProcess.options.ignore_cols_all_zeros.Value   =  0;
    
    sProcess.options.ignore_rows_all_zeros.Comment = 'Ignore all-zero rows';
    sProcess.options.ignore_rows_all_zeros.Type    = 'checkbox';
    sProcess.options.ignore_rows_all_zeros.Value   =  0;
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
if  sProcess.options.ignore_rows_all_zeros.Value
    rows_to_keep = any(DataMat.Value ~= 0, 2)';
else
    rows_to_keep = 1:size(DataMat.Value, 1);
end
if  sProcess.options.ignore_cols_all_zeros.Value
    cols_to_keep = any(DataMat.Value ~= 0, 1);
else
    cols_to_keep = 1:size(DataMat.Value, 2);
end

[row_names, col_names] = process_nst_concat_matrices('get_axes_info', DataMat);
table_out = array2table(DataMat.Value(rows_to_keep, cols_to_keep), ...
                        'VariableNames', DataMat.ColNames(cols_to_keep), ...
                        'RowNames', DataMat.RowNames(rows_to_keep));
writetable(table_out, sProcess.options.csv_file.Value{1}, ...
           'Delimiter', ',', 'WriteRowNames', 1);
end
