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


sProcess.options.prefixes.Comment = 'Row or Column prefixes (comma-separated)';
sProcess.options.prefixes.Type = 'text';
sProcess.options.prefixes.Value = '';

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

stacking_types = get_stacking_types();
stacking_type_int = sProcess.options.stacking_type.Value{1};
stacking_type_str = sProcess.options.stacking_type.Value{2}{stacking_type_int};

[varying_comment, common_prefix, common_suffix] = str_remove_common({sInputs.Comment}, 1);

if ~isempty(sProcess.options.prefixes.Value)
    prefixes = cellfun(@strtrim, strsplit(sProcess.options.prefixes.Value, ','),...
        'UniformOutput', 0);
    
    if length(sInputs) ~= length(prefixes)
        error('Number of prefixes not equal to number of inputs');
    end
elseif ~isempty(varying_comment) && length(varying_comment)==length(sInputs)
    prefixes = cellfun(@(c) [c '_'], varying_comment, 'UniformOutput', 0);
else
    prefixes = repmat({''}, 1, length(sInputs));
end

switch stacking_type_int
    case stacking_types.row
    case stacking_types.column
    otherwise
        error('Invalid stacking type');
end

% Load the first file, as the reference
MatRef = in_bst_matrix(sInputs(1).FileName);
[row_names_ref, col_names_ref, sfreq_ref, time_ref] = get_axes_info(MatRef);
nb_rows_ref = length(row_names_ref);
nb_cols_ref = length(col_names_ref);

% Allocate concatenated items
switch stacking_type_int
    case stacking_types.row
        nb_cats = nb_rows_ref * length(sInputs);
        values = zeros(nb_cats, nb_cols_ref);
        row_names = cell(1, nb_cats);
        row_names(1:nb_rows_ref) = cellfun(@(c) [prefixes{1} c], row_names_ref, 'UniformOutput', 0);
        col_names = col_names_ref;
        time = zeros(1, nb_cats);
        time(1:nb_rows_ref) = time_ref;
    case stacking_types.column
        nb_cats = nb_cols_ref * length(sInputs);
        values = zeros(nb_rows_ref, nb_cats);
        row_names = row_names_ref;
        col_names = cell(1, nb_cats);
        col_names(1:nb_cols_ref) =  cellfun(@(c) [prefixes{1} c], col_names_ref, 'UniformOutput', 0);
        time = time_ref;
end
values_std = zeros(size(values));

values(1:size(MatRef.Value, 1), 1:size(MatRef.Value, 2)) = MatRef.Value;
values_std(1:size(MatRef.Std, 1), 1:size(MatRef.Std, 2)) = MatRef.Std;

MatNew = db_template('matrix');

% Set history field
MatNew = bst_history('add', MatNew, 'concat', ...
                     sprintf('Contatenate %ss from files:', stacking_type_str));
MatNew = bst_history('add', MatNew, 'concat', [' - ' sInputs(1).FileName]);

% Loop over remaining input files
for iInput = 2:length(sInputs)
    % Load the next file
    MatToCat = in_bst_matrix(sInputs(iInput).FileName);
    
    % Check consistency
    [row_names_new, col_names_new, sfreq_new, time_new] = get_axes_info(MatToCat);
    
    if stacking_type_int==stacking_types.row && ...
            (length(col_names_new) ~= length(col_names_ref) || ...
            ~all(strcmp(col_names_new, col_names_new)))
        bst_report('Error', sProcess, sInputs(iInput), ['This file has columns incompatible with the first one: "' sInputs(iInput).FileName '".']);
        return;
    elseif stacking_type_int==stacking_types.column && ...
            (length(row_names_new) ~= length(row_names_ref) || ...
            ~all(strcmp(row_names_new, row_names_new)))
        bst_report('Error', sProcess, sInputs(iInput), ['This file has rows incompatible with the first one: "' sInputs(iInput).FileName '".']);
        return;
    end
    
    % Concatenate the events
    if isfield(MatToCat, 'Events') && ~isempty(MatToCat.Events)
        if isempty(sfreq_new)
            bst_report('Warning', sProcess, sInputs(iInput), ...
                ['Events of this file are skipped because time is not defined: "' sInputs(iInput).FileName '".']);
        else
            % Add the events to the new file
            if isempty(MatRef.Events)
                MatRef.Events = MatToCat.Events;
            else
                % Trick import_events() to work for event concatenation
                sFile.events = MatNew.Events;
                sFile.prop.sfreq = sfreq_new;
                sFile = import_events(sFile, [], MatToCat.Events);
                MatNew.Events = sFile.events;
            end
        end
    end
    
    % Concatenate the data matrices
    switch stacking_type_int
        case stacking_types.row
            irow = (nb_rows_ref*(iInput-1)+1):nb_rows_ref*iInput;
            values(irow, :) = MatToCat.Value;
            if isfield(MatToCat, 'Std') && ~isempty(MatToCat.Std)
                values_std(irow, :) = MatToCat.Std;
            end
            row_names(irow) = cellfun(@(c) [prefixes{iInput} c], row_names_new, 'UniformOutput', 0);
            
            if stacking_type == stacking_types.row && ~isempty(time_ref)
                % Concatenate Time axis 
                time(irow) = time_new + time(irow(end)-1);
            end
            
        case stacking_types.column
            icol = (nb_cols_ref*(iInput-1)+1):nb_cols_ref*iInput;
            values(:, icol) = MatToCat.Value;
            col_names(icol) =  cellfun(@(c) [prefixes{iInput} c], col_names_new, 'UniformOutput', 0);
            if isfield(MatToCat, 'Std') && ~isempty(MatToCat.Std)
                values_std(:, icol) = MatToCat.Std;
            end
    end
    
    % History field
    MatNew = bst_history('add', MatNew, 'concat', [' - ' sInputs(iInput).FileName]);
end

MatNew.Value = values;
MatNew.Std = values_std;
MatNew.RowNames = row_names;
MatNew.Time = time;
MatNew.ColNames = col_names;
MatNew.Description = col_names;
MatNew.DisplayUnits = MatRef.DisplayUnits; % TODO check unit consistency across inputs

% Output file tag
fileTag = 'matrix_concat';

% Set comment
MatNew.Comment = [common_prefix, common_suffix, ...
                  sprintf(' | concat (%d) %ss', nb_cats, stacking_type_str)];
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

end
