function varargout = process_nst_concat_matrices( varargin )
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
sProcess.Index       = 1400;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'matrix'};
sProcess.OutputTypes = {'matrix'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 2;

% Definition of the options
sProcess.options.stacking_type.Comment = 'Stack by:';
sProcess.options.stacking_type.Type    = 'combobox';
choices = get_stacking_types();
sProcess.options.stacking_type.Value   = {choices.column,...
    fieldnames(choices)};


% TODO: use dedicated process to add a prefix to col or row names
% -> cannot easily manage it here since nb of rows/cols can be variable
% sProcess.options.prefixes.Comment = 'Row or Column prefixes (comma-separated)';
% sProcess.options.prefixes.Type = 'text';
% sProcess.options.prefixes.Value = '';

end

function stypes = get_stacking_types()
stypes.row = 1;
stypes.column = 2;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};

stacking_types = get_stacking_types();
stacking_type_int = sProcess.options.stacking_type.Value{1};
stacking_type_str = sProcess.options.stacking_type.Value{2}{stacking_type_int};

[varying_comment, common_prefix, common_suffix] = str_remove_common({sInputs.Comment}, 1);
varying_comment(cellfun(@isempty, varying_comment)) = {''};

% if ~isempty(sProcess.options.prefixes.Value)
%     prefixes = cellfun(@strtrim, strsplit(sProcess.options.prefixes.Value, ','),...
%         'UniformOutput', 0);
%     
%     if length(sInputs) ~= length(prefixes)
%         error('Number of prefixes not equal to number of inputs');
%     end
% elseif ~isempty(varying_comment) && length(varying_comment)==length(sInputs)
%     prefixes = cellfun(@(c) [c '_'], varying_comment, 'UniformOutput', 0);
% else
%     prefixes = repmat({''}, 1, length(sInputs));
% end

switch stacking_type_int
    case stacking_types.row
    case stacking_types.column
    otherwise
        error('Invalid stacking type');
end

% Load the first file, as the reference
MatRef = in_bst_matrix(sInputs(1).FileName);
[row_names_ref, col_names_ref] = get_axes_info(MatRef);
row_names = row_names_ref;
col_names = col_names_ref;

% if ~isempty(MatRef.Time)
%    warning(sprintf('Time field not empty in "%s". It will be dropped.\n Consider using field RowName.', ...
%                    sInputs(1).FileName));
% end
% 
% if ~isempty(MatRef.Description)
%    warning(sprintf('Description field not empty in "%s". It will be dropped.\n Consider using field ColName.', ...
%                    sInputs(1).FileName));
% end

if isfield(MatRef, 'Events') && ~isempty(MatRef.Events)
   warning('Events field not empty in "%s". It will be dropped.', sInputs(1).FileName);
end
    
values = MatRef.Value;
values_std = MatRef.Std;
time = MatRef.Time;
description = MatRef.Description;

MatNew = db_template('matrix');

% Set history field
MatNew = bst_history('add', MatNew, 'concat', ...
                     sprintf('Contatenate %ss from files:', stacking_type_str));
MatNew = bst_history('add', MatNew, 'concat', [' - ' sInputs(1).FileName]);

for iInput=2:length(sInputs)
    MatToCat = in_bst_matrix(sInputs(iInput).FileName);

%     if ~isempty(MatToCat.Time)
%        warning(sprintf('Time field not empty in "%s". It will be dropped.\n Consider using field RowName.', ...
%                        sInputs(iInput).FileName));
%     end
% 
%     if ~isempty(MatToCat.Description)
%        warning(sprintf('Description field not empty in "%s". It will be dropped.\n Consider using field ColName.', ...
%                        sInputs(iInput).FileName));
%     end    
    
    if isfield(MatToCat, 'Events') && ~isempty(MatToCat.Events)
       warning(sprintf('Events field not empty in "%s". It will be dropped.', ...
                       sInputs(iInput).FileName));
    end
    
    % Check consistency
    [row_names_new, col_names_new] = get_axes_info(MatToCat);

    if stacking_type_int==stacking_types.row
        if (length(col_names_new) ~= length(col_names_ref) || ...
            ~all(strcmp(col_names_new, col_names_ref)))
            msg = ['This file has columns incompatible with the first one: "' sInputs(iInput).FileName '".'];
            throw(MException('Nirstorm:IncompatibleMatrices', msg)); 
        end
        
%         common_rows = ismember(row_names_new, row_names_ref);
%         if any(common_rows)
%             bst_report('Error', sProcess, sInputs(iInput), ['This file has common rows with the first one: "' sInputs(iInput).FileName '".']);
%             return;
%         end
        values = [values ; MatToCat.Value];
        values_std = [values_std ; MatToCat.Std];
        row_names = [row_names row_names_new];
        time = [time MatToCat.Time];
    elseif stacking_type_int==stacking_types.column
        if (length(row_names_new) ~= length(row_names_ref) || ...
                ~all(strcmp(row_names_new, row_names_ref)))
            msg = ['This file has rows incompatible with the first one: "' sInputs(iInput).FileName '".'];
            throw(MException('Nirstorm:IncompatibleMatrices', msg));
        end
%         common_cols = ismember(col_names_new, col_names_ref);
%         if any(common_cols)
%             bst_report('Error', sProcess, sInputs(iInput), ['This file has common columns with the first one: "' sInputs(iInput).FileName '".']);
%             return;
%         end
        
        values = [values MatToCat.Value];
        values_std = [values_std MatToCat.Std];
        col_names = [col_names col_names_new];
        description = [description MatToCat.Description];
    end
       
    MatNew = bst_history('add', MatNew, 'concat', [' - ' sInputs(iInput).FileName]);
    
end

MatNew.Value = values;
MatNew.Std = values_std;
MatNew.RowNames = row_names;
MatNew.Time = time;
MatNew.ColNames = col_names;
MatNew.Description = description;
MatNew.DisplayUnits = MatRef.DisplayUnits; % TODO check unit consistency across inputs

% Output file tag
fileTag = 'matrix_concat';

% Set comment
MatNew.Comment = [common_prefix, strjoin(varying_comment, ' '), common_suffix, ...
                  sprintf(' | concat %ss (%d files) ', stacking_type_str, length(sInputs))];
% Get output filename
OutputFiles{1} = bst_process('GetNewFilename', bst_fileparts(sInputs(1).FileName), fileTag);
% Save file
bst_save(OutputFiles{1}, MatNew, 'v6');
% Register in database
db_add_data(sInputs(1).iStudy, OutputFiles{1}, MatNew);
end

function [row_names, col_names, sfreq, time] = get_axes_info(DataMat)

if ~isempty(DataMat.Time) && isnumeric(DataMat.Time)
    time = DataMat.Time;
    sfreq = 1 ./ (DataMat.Time(2) - DataMat.Time(1));
else
    time = []; 
    sfreq = [];
end

if isfield(DataMat, 'RowNames') && ~isempty(DataMat.RowNames)
    row_names = DataMat.RowNames;
elseif ~isempty(DataMat.Time) && isnumeric(DataMat.Time)
    row_names = arrayfun(@(n) num2str(n), DataMat.Time, 'UniformOutput', 0);
else 
    error('Cannot resolve matrix row names');
end

if isfield(DataMat, 'ColNames') && ~isempty(DataMat.ColNames)
    col_names = DataMat.ColNames;
elseif ~isempty(DataMat.Description)
    col_names = DataMat.Description;
    if isnumeric(col_names)
        col_names = arrayfun(@(n) num2str(n), col_names, 'UniformOutput', 0);
    end 
else 
    error('Cannot resolve matrix column names');
end

row_names = line_vector(row_names);
col_names = line_vector(col_names);

end


function v = line_vector(v)
if size(v, 2) == 1 && size(v, 1) > 1
    v = v';
end
end