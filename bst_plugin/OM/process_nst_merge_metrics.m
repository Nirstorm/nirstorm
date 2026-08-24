function varargout = process_nst_merge_metrics( varargin )
% PROCESS_NST_MERGE_METRICS: Fusionne plusieurs tableaux de métriques
% et sauvegarde le tableau unique résultant dans un dossier 'metrics' du sujet.
% Utilise la convention officielle 2D Description de Brainstorm (Commit c3b18b5).
%
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2026 University of Southern California & McGill University
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
% Authors: Edouard Delaire, Jean-Eudes Bornert 2026

eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description du process
    sProcess.Comment     = 'Merge montage metrics';
    sProcess.FileTag     = '';
    
    sProcess.Category    = 'Custom'; 
    sProcess.SubGroup    = {'NIRS', 'Optimal Montage'};
    sProcess.Index       = 1109;
    
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'matrix'}; 
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
    sProcess.isSeparator = 0;
    
    sProcess.options.output_name.Comment = 'Output table name:';
    sProcess.options.output_name.Type    = 'text';
    sProcess.options.output_name.Value   = 'Merged_Metrics_Table';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

%% ===== RUN ======
function OutputFiles = Run(sProcess, sInput)
    OutputFiles = {};
    
    TargetCondition = 'metrics';
    SubjectName = sInput(1).SubjectName;
    
    all_values = [];
    all_rownames = {};
    col_names = {};
    display_units = '';
    
    bst_progress('start', 'Merging metrics', 'Reading and merging tables...', 0, length(sInput));
    
    for iFile = 1:length(sInput)
        
        sMat = in_bst(sInput(iFile).FileName);
        
        val = sMat.Value;
        [nRows, nCols] = size(val);
        desc = sMat.Description;
        
        current_cols = cell(1, nCols);
        current_rows = cell(nRows, 1);
        
        if ~isempty(desc) && iscell(desc) && size(desc,1) >= 1 && size(desc,2) >= 1
            first_elem = desc{1,1};
            if iscell(first_elem) && length(first_elem) == 2
                current_rows{1} = first_elem{1};
                current_cols{1} = first_elem{2};
            else
                current_rows{1} = sInput(iFile).Condition;
                if ischar(first_elem) || isstring(first_elem)
                    current_cols{1} = first_elem;
                else
                    current_cols{1} = 'Col_1';
                end
            end
            
            if nCols > 1 && size(desc,2) >= nCols
                current_cols(2:nCols) = desc(1, 2:nCols);
            end
            
            if nRows > 1 && size(desc,1) >= nRows
                current_rows(2:nRows) = desc(2:nRows, 1);
            end
        else
            current_cols = arrayfun(@(n) sprintf('Col_%d', n), 1:nCols, 'UniformOutput', 0);
            current_rows = arrayfun(@(n) sprintf('%s_%d', sInput(iFile).Condition, n), 1:nRows, 'UniformOutput', 0);
        end
        
        current_cols = line_vector(current_cols);
        
        if iFile == 1
            col_names = current_cols;
            if isfield(sMat, 'DisplayUnits')
                display_units = sMat.DisplayUnits;
            end
        else
            if ~isequal(col_names, current_cols)
                bst_error(['Inconsistencies in the column headers between the selected tables : ' 10 ...
                           'It is not possible to merge tables with different structures.']);
                return;
            end
        end
        
        if size(val, 1) ~= length(current_rows)
            val = reshape(val, length(current_rows), length(col_names));
        end
        
        all_values = [all_values; val];
        all_rownames = [all_rownames; current_rows(:)];
        
        bst_progress('set', iFile);
    end
    
    bst_progress('stop');
    
    all_rownames = matlab.lang.makeUniqueStrings(all_rownames);
    
    merged_table = array2table(all_values, 'VariableNames', col_names);
    merged_table.Properties.RowNames = all_rownames;
    
    if isfield(sProcess.options, 'output_name') && ~isempty(sProcess.options.output_name.Value)
        CommentName = sProcess.options.output_name.Value;
    else
        CommentName = 'Merged Metrics Table';
    end

    sFile = save_table(merged_table, SubjectName, TargetCondition, CommentName, struct(), display_units);
    
    OutputFiles = {sFile}; 
end

%% ===== HELPER FUNCTIONS =====
function sFile = save_table(t, subject_name, condition, comment, extra, displayUnits)

    % Save a table as a Brainstorm matrix item using the official 2D Description format (Commit c3b18b5)
    if nargin < 5 
        extra = struct();
    end
    if nargin < 6
        displayUnits = '';
    end
    
    MatNew = db_template('matrixmat');
    MatNew.Value = table2array(t);
    [nRows, nCols] = size(MatNew.Value);
    
    MatNew.Std = [];
    MatNew.Time = []; 
    MatNew.DisplayUnits = displayUnits;
    MatNew.Comment = comment;
    
    if ~isempty(t.Properties.RowNames)
        rowNames = t.Properties.RowNames(:);
    else
        rowNames = arrayfun(@(n) sprintf('row_%d', n), 1:nRows, 'UniformOutput', 0)';
    end
    colNames = t.Properties.VariableNames(:)';
    
    MatNew.Description = cell(nRows, nCols);
    
    MatNew.Description{1, 1} = {rowNames{1}, colNames{1}};
    
    if nCols > 1
        MatNew.Description(1, 2:nCols) = colNames(2:nCols);
    end
    
    if nRows > 1
        MatNew.Description(2:nRows, 1) = rowNames(2:nRows);
    end
    
    sSubject = bst_get('Subject', subject_name, 1);
    if isempty(sSubject)
        db_add_subject(subject_name, []);
    end
    
    [sStudy, iStudy] = bst_get('StudyWithCondition', [subject_name '/' condition]);
    if isempty(sStudy) 
        iStudy = db_add_condition(subject_name, condition);
        sStudy = bst_get('Study', iStudy);
    end
    
    extra_fields = fieldnames(extra);
    for ifield = 1:length(extra_fields)
        MatNew.(extra_fields{ifield}) = extra.(extra_fields{ifield});
    end
    
    sFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_table');
    bst_save(sFile, MatNew, 'v6');
    
    db_add_data(iStudy, sFile, MatNew);
end

function v = line_vector(v)
    if size(v, 2) == 1 && size(v, 1) > 1
        v = v';
    end
end