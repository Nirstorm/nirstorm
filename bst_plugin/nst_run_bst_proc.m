function sFilesOut = nst_run_bst_proc(out_items_names, force_redo, ProcessName, sFiles, sFiles2, varargin)
% Allow to run a brainstorm process only if expected outputs do not already exist 
% or if recomputation is forced.
% Outputs have predefined names as given by out_items_names. Each output
% will be unique. This overrides the default behaviour of brainstorm 
% to always run the process and add a suffix to the output name if
% an item with the same name already exists.
% 
% nst_run_bst_proc forces the called process to have unique outputs. If outputs
% already exist, the process will not be executed.
% Unless force_redo=1. In this case, the exisiting outputs will be deleted
% prior to executing the process.
% 
% Note that the definition of the output names (Comment field) is
% overriden and no longer determined by the called process. Therefore this
% function is rather suited only for processes with a fixed predictible number 
% of outputs.
%
% WARNING: this is a helper for "simple" processes producing predictible
% outputs in the same condition folder as the input data.
% 
% WARNING: works only for functional data
%
% Input:
%    - out_items_names (str or cell array of str):
%          Output item comment(s). Must be consistent with the number of outputs 
%          returned by process "ProcessName".
%          The output condition can be specified as: '<condition>/<comment>'
%          The slash character is not allowed in <condition> and <comment>.
%          If no condition is given, the one from the given sFiles is taken.
%    - force_redo (bool):
%          Flag to force recomputation.
%    - sFiles (cell array of file names or array of bst process input structures):
%          Input data as in the 3rd arg of bst_process('CallProcess', ProcessName, sFiles, sFiles2).
%    - sFiles2 (cell array of file names or array of bst process input structures):
%          Input data as in the 4th arg of bst_process('CallProcess', ProcessName, sFiles, sFiles2). 
%    - varagin: 
%          options passed to the brainstorm process.
%
%% Example
% Basic brainstorm usage is:
%    sFilesOut = bst_process('CallProcess', 'process', sFiles, sFiles2, ...
%                            'proc_option1', proc_option1_val,...
%                            'proc_option2', proc_option2_val);
% This creates items whose names are determined by 'process', for instance
% "result" and "result_other".
% If 'process' is called again with the same parameters, new outputs will
% be created: "result_02" and "result_other_02".
%
% To avoid this, the example above can be translated to:
%    sFilesOut = nst_run_bst_proc({'result', 'result_other'}, 0, ...
%                                 'process', sFiles, sFiles2, ...
%                                 'proc_option1', proc_option1_val,...
%                                 'proc_option2', proc_option2_val);
% In this case if "result" and "result_other" already exist,
% the process will not be executed. If they don't exist, the process is
% executed and its outputs are renamed to "result" and "result_other".
% 
% %% Force recomputation %%
%    redo = 1;
%    sFilesOut = nst_run_bst_proc({'result', 'result_other'}, redo, ...
%                                 'process', sFiles, sFiles2, ...
%                                 'proc_option1', proc_option1_val,...
%                                 'proc_option2', proc_option2_val);
% In this case, if "result" and "result_other" already exist, they will be deleted
% before running the process.
%
% %% Specify output condition folder %%
%    sFilesOut = nst_run_bst_proc({'new_cond/result', 'new_cond/result_other'}, redo, ...
%                                 'process', sFiles, sFiles2, ...
%                                 'proc_option1', proc_option1_val,...
%                                 'proc_option2', proc_option2_val);
% 
% Outputs 'result' and 'result_other' will be placed in folder 'new_cond'
% (created if needed).

out_name_re = '(?<condition>.*/)?(?<comment>[^/]+)';

%% Check inputs
if nargin < 5
    sFiles2 = [];
end

if ischar(out_items_names)
    out_items_names = {out_items_names};
end
if ~iscell(out_items_names) || ~all(cellfun(@ischar, out_items_names))
    throw(MException('Nirstorm:BadArgType', 'out_items_names must be str or cell array of str.')); 
end

if ~iscell(out_items_names) || any(cellfun(@isempty, regexp(out_items_names, out_name_re, 'match')))
    throw(MException('Nirstorm:BadArgType', ...
                     'output item name must be formated as "output_name" or "condition/output_name".')); 
end

if ~(force_redo == 1 || force_redo == 0)
   throw(MException('Nirstorm:BadArgType', 'force_redo must be 0 or 1')); 
end

if ~ischar(ProcessName)
   throw(MException('Nirstorm:BadArgType', 'ProcessName must be str'));  
end

if ~isstruct(sFiles)
    sInputs = bst_process('GetInputStruct', sFiles);
else
    sInputs = sFiles;
end
assert(length(unique({sInputs.SubjectName}))==1);
subject_name = sInputs(1).SubjectName;

%% Parse output definitions (extract condition if given)
outputs = struct();
for i_item=1:length(out_items_names)
    toks = regexp(out_items_names{i_item}, out_name_re, 'names');
    assert(length(toks)==1);
    outputs(i_item).condition = replace(toks.condition, '/', '');
    outputs(i_item).comment = toks.comment;
    if isempty(outputs(i_item).condition)
        outputs(i_item).condition =  sInputs(1).Condition;
    end
end

%% Look for existing outputs
sFilesOut = {};
duplicates = {};

for i_item=1:length(outputs)
    % Manage target condition folder
    dest_condition_folder = bst_fullfile(subject_name, outputs(i_item).condition);
    [sStudy, iStudy] = bst_get('StudyWithCondition', dest_condition_folder);
    if isempty(sStudy) % condition does not exist -> create it
        iStudy = db_add_condition(subject_name, outputs(i_item).condition);
        sStudy = bst_get('Study', iStudy);
    end
    outputs(i_item).sStudy = sStudy;
    outputs(1).iStudy = iStudy;
    
    % Check if output exists in specified condition
    [selected_files, file_type] = nst_get_bst_func_files(subject_name, outputs(i_item).condition, outputs(i_item).comment);
    if ~isempty(selected_files) && ~ischar(selected_files) && length(selected_files) > 1
        duplicates{end+1} = [outputs(i_item).condition '/' outputs(i_item).comment];
    end
    if ~isempty(file_type)
        sFilesOut_types{i_item} = file_type;
    else
        sFilesOut_types{i_item} = '';
    end
    sFilesOut{i_item} = selected_files;
end
if ~isempty(duplicates)
    bst_error(sprintf('Cannot safely manage unique outputs. Found duplicate items: %s', strjoin(duplicates, ', ')));
    sFilesOut = {};
    return;
end
existing = cellfun(@(s) ~isempty(s), sFilesOut);

%% Run the process if needed
if any(~existing) || force_redo
    if any(existing)
        if strcmp(sFilesOut_types{1}, 'HeadModel')
            assert(length(sFilesOut_types) == 1);
            prev_iHeadModel = strcmp({outputs(1).sStudy.HeadModel.FileName}, sFilesOut{1});
            outputs(1).sStudy = delete_head_model(outputs(1).sStudy, outputs(1).iStudy, prev_iHeadModel);
        else
            bst_process('CallProcess', 'process_delete', sFilesOut, [], ...
                       'target', 1);
        end
        bst_report('Info', ProcessName, sFiles, ...
                   sprintf('Force redo - removed previous result(s): %s', strjoin(sFilesOut, ', ')) );
    end
    
    % Special case for head model which is not returned in sFilesOut
    % -> keep track of iHeadmodel 
    prev_iHeadmodel = sStudy.iHeadModel;
    
    % Call the process
    sFilesOut = bst_process('CallProcess', ProcessName, sFiles, sFiles2, varargin{:});
    
    % Check if process created a new head model
    sStudy = bst_get('Study', sInputs(1).iStudy);
    new_iHeadModel = setdiff(sStudy.iHeadModel, prev_iHeadmodel);
    assert(length(new_iHeadModel) <= 1); %just a safe-guard, should always be the case
    
    if isstruct(sFilesOut)
        sFilesOut = {sFilesOut.FileName};
    end
    
    % Check outputs consistency and rename them
    if isempty(sFilesOut) && ~isempty(new_iHeadModel) % special case for head model
        if length(outputs) ~= 1
            delete_head_model(sStudy, sInputs(1).iStudy, new_iHeadModel);
            bst_error(sprintf('Expected %d outputs but process produced only one head model.\n', ...
                              length(out_items_names)));
            sFilesOut = {};
            return;
        end
        if sInputs(1).iStudy == outputs(1).iStudy
            rename_head_model(sStudy, sInputs(1).iStudy, new_iHeadModel, out_items_names{1});
        else
            error('Moving of head model to new condition not supported');
        end
    else
        if length(sFilesOut) ~= length(outputs) && ~strcmp(sFilesOut_types{1}, 'HeadModel')
            bst_process('CallProcess', 'process_delete', sFilesOut, [], ...
                        'target', 1); 
            bst_error(sprintf('Expected %d outputs but process produced %d.\n', ...
                              length(outputs), length(sFilesOut)));
            sFilesOut = {};
            return;
        end
        for i_item=1:length(outputs)
            sOutRenamed = bst_process('CallProcess', 'process_set_comment', sFilesOut{i_item}, [], ...
                               'tag', outputs(i_item).comment, ...
                               'isindex', 0);

            if ~strcmp(sOutRenamed.Condition, outputs(i_item).condition)
                sOut = bst_process('CallProcess', 'process_movefile', sOutRenamed, [], ...
                                   'subjectname', subject_name, ...
                                   'folder', outputs(i_item).condition);
                [sChan, iChan] = bst_get('ChannelForStudy',   outputs(i_item).iStudy);
                
                % Copy channel file along with data if needed
                if ismember(sOut.FileType, {'data', 'pdata'}) && isempty(sChan)
                    ChannelMat = in_bst_channel(sOutRenamed.ChannelFile);
                    db_set_channel(iChan, ChannelMat, 0, 0);
                end
            else
                sOut = sOutRenamed;
            end               
            if isstruct(sOut)
                sFilesOut{i_item} = sOut.FileName;
            else
                sFilesOut{i_item} = sOut;
            end
        end
    end
else % no need to run the process
    bst_report('Info', ProcessName, sFiles, ...
               sprintf('Skipped execution of %s. Outputs found.', ProcessName));
    if strcmp(sFilesOut_types{1}, 'HeadModel')
        assert(length(sFilesOut) == 1);
        sFilesOut = {}; % Do not return any output for head model computation to comply with brainstorm way
    end
end

if length(sFilesOut) == 1
    sFilesOut = sFilesOut{1};
end

end

function sStudy = delete_head_model(sStudy, iStudy, iHeadModelDel)

% From node_delete.m / case 'headmodel'

% Remove files descriptions from database
sStudy.HeadModel(iHeadModelDel) = [];
% Update default headmodel
nbHeadModel = length(sStudy.HeadModel);
if (nbHeadModel <= 0)
    sStudy.iHeadModel = [];
elseif (nbHeadModel == 1)
    sStudy.iHeadModel = 1;
elseif (sStudy.iHeadModel > nbHeadModel)
    sStudy.iHeadModel = nbHeadModel;
else
    % Do not change iHeadModel
end
% Study was modified
bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);
db_save();
end

function sStudy = rename_head_model(sStudy, iStudy, iHeadModelRename, new_name)
sStudy.HeadModel(iHeadModelRename).Comment = new_name;
bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);
db_save();
end