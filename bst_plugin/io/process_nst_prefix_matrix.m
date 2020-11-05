function varargout = process_nst_prefix_matrix( varargin )
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2018 University of Southern California & McGill University
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
% Authors: Thomas Vincent, 2019

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Concatenate matrices';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS - wip';
sProcess.Index       = 1904;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'matrix'};
sProcess.OutputTypes = {'matrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

% Definition of the options


sProcess.options.row_prefixes.Comment = 'Row prefix(es) (comma-separated)';
sProcess.options.row_prefixes.Type = 'text';
sProcess.options.row_prefixes.Value = '';

sProcess.options.col_prefixes.Comment = 'Column prefix(es) (comma-separated)';
sProcess.options.col_prefixes.Type = 'text';
sProcess.options.col_prefixes.Value = '';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};

row_prefixes = cellfun(@strtrim, strsplit(sProcess.options.row_prefixes.Value, ','),...
                       'UniformOutput', 0);

col_prefixes = cellfun(@strtrim, strsplit(sProcess.options.col_prefixes.Value, ','),...
                       'UniformOutput', 0);                       
if length(row_prefixes) == 1 && isempty(row_prefixes{1}) && ...
        length(col_prefixes) == 1 && isempty(col_prefixes{1})
    return
end
                                  
MatData = in_bst_matrix(sInputs(1).FileName);
if ~isfield(MatData, 'RowNames')
    error('Field RowNames not found in given matrix');
end

if ~isfield(MatData, 'ColNames')
    error('Field ColNames not found in given matrix');
end

if size(MatData.Value, 1) ~= length(row_prefixes)
    if length(row_prefixes) > 1
        error('Number of row prefixes (%d) not equal to number of rows (%d)', ...
              length(row_prefixes), size(MatData.Value, 1));
    else
        row_prefixes = repmat(row_prefixes, 1, size(MatData.Value, 1));
    end
end

if size(MatData.Value, 2) ~= length(col_prefixes)
    if length(col_prefixes) > 1
        error('Number of column prefixes (%d) not equal to number of columns (%d)', ...
              length(row_prefixes), size(MatData.Value, 2));
    else
        col_prefixes = repmat(col_prefixes, 1, size(MatData.Value, 2));
    end
end

% Set history field
if length(col_prefixes) > 1 || ~isempty(col_prefixes{1})
    MatData = bst_history('add', MatData, 'Renamed columns');
end

if length(row_prefixes) > 1 || ~isempty(row_prefixes{1})
    MatData = bst_history('add', MatData, 'Renamed rows');
end

for irow=1:length(MatData.RowNames)
    MatData.RowNames{irow} = [row_prefixes{irow} MatData.RowNames{irow}];
end
for icol=1:length(MatData.ColNames)
    MatData.ColNames{icol} = [col_prefixes{icol} MatData.ColNames{icol}];
end

MatData.Description = MatData.ColNames;

%TODO: allow saving in new file

% Get output filename
OutputFiles{1} = file_fullpath(sInputs(1).FileName);
% Save file
bst_save(OutputFiles{1}, MatData, 'v6');
end
