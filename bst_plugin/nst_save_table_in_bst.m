function sFile = nst_save_table_in_bst(t, subject_name, condition, comment, extra, displayUnits)
% Save a table as a matrix brainstorm item:
%   Field Std is not used
%   Fields Time and Description are not used because their interpretation by
%   brainstorm is not straightfoward.
%   Description sometimes refer to row names and sometimes to column names
%   depending on the matrix dimension and if Time is available.
%   Here new fields RowNames and ColNames are added to store them.
%
%

if nargin < 5 
    extra = struct();
end

if nargin < 6
    displayUnits = '';
end

MatNew = db_template('matrix');
MatNew.Value = table2array(t);
nrows = size(MatNew.Value, 1);
MatNew.Std = [];
MatNew.Time = [];
MatNew.Description = {};
MatNew.DisplayUnits = displayUnits;
MatNew.Comment = comment;

if ~isempty(t.Row)
    MatNew.RowNames = t.Row;
else
    MatNew.RowNames = arrayfun(@(n) sprintf('row_%d', n), 1:nrows, 'UniformOutput', 0);
end
MatNew.RowNames = line_vector(MatNew.RowNames);
MatNew.ColNames = line_vector(t.Properties.VariableNames);


sSubject = bst_get('Subject', subject_name, 1);
if isempty(sSubject)
    db_add_subject(subject_name, []);
end
[sStudy, iStudy] = bst_get('StudyWithCondition', [subject_name '/' condition]);
if isempty(sStudy) 
    iStudy = db_add_condition(subject_name, condition);
    sStudy = bst_get('Study', iStudy);
end

% Add extra fields
 extra_fields = fieldnames(extra);
 for ifield = 1:length(extra_fields)
%      assert(~isfield(MatNew, extra_fields{ifield}));
     MatNew.(extra_fields{ifield}) = extra.(extra_fields{ifield});
 end

sFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'matrix_table');
% Save file
bst_save(sFile, MatNew, 'v6');
% Register in database
db_add_data(iStudy, sFile, MatNew);
end

function v = line_vector(v)
if size(v, 2) == 1 && size(v, 1) > 1
    v = v';
end
end