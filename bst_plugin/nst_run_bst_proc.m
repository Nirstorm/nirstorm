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
%          Output item(s) names(s), looked for in the condition folder of sFile. 
%          Must have the same length as the actual number of outputs 
%          given by process "ProcessName".
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
% Calling:
%    sFilesOut = bst_process('CallProcess', 'process', sFiles, sFiles2, ...
%                            'proc_option1', proc_option1_val,...
%                            'proc_option2', proc_option2_val);
% will create items whose names are determined by the process, for example
% "result" and "result_other".
% If the process is called again with the same parameters, new outputs will
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
% If ones wants to force recomputation:
%    sFilesOut = nst_run_bst_proc({'result', 'result_other'}, 0, ...
%                                 'process', sFiles, sFiles2, ...
%                                 'proc_option1', proc_option1_val,...
%                                 'proc_option2', proc_option2_val);
% In this case, if "result" and "result_other" already exist, they will be deleted
% and the process will be run.
%
% TODO: handle anatomy outputs

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

%% Look for existing outputs
sFilesOut = {};
duplicates = {};
for i_item=1:length(out_items_names)
    selected_files = nst_get_bst_func_files(sInputs(1).SubjectName, sInputs(1).Condition, out_items_names{i_item});
    if ~isempty(selected_files) && ~ischar(selected_files) && length(selected_files) > 1
        duplicates{end+1} = out_items_names{i_item};
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
        bst_report('Info', ProcessName, sFiles, ...
                   sprintf('Force redo - removing previous result(s): %s', strjoin(sFilesOut, ', ')) );

        bst_process('CallProcess', 'process_delete', sFilesOut, [], ...
                    'target', 1); 
    end
    sFilesOut = bst_process('CallProcess', ProcessName, sFiles, sFiles2, varargin{:});
    if isstruct(sFilesOut)
        sFilesOut = {sFilesOut.FileName};
    end
    
    if length(sFilesOut) ~= length(out_items_names)
        bst_process('CallProcess', 'process_delete', sFilesOut, [], ...
                    'target', 1); 
        bst_error(sprintf('Expected %d outputs but process produced %d.\n', ...
                          length(out_items_names), length(sFilesOut)));
        sFilesOut = {};
        return;
    end
    for i_item=1:length(out_items_names)
        sOut = bst_process('CallProcess', 'process_set_comment', sFilesOut{i_item}, [], ...
                                        'tag', out_items_names{i_item}, ...
                                        'isindex', 0);
        if isstruct(sOut)
            sFilesOut{i_item} = sOut.FileName;
        else
            sFilesOut{i_item} = sOut;
        end
    end
else
    bst_report('Info', ProcessName, sFiles, ...
               sprintf('Skipped execution of %s. Outputs found.', ProcessName));
end

if length(sFilesOut) == 1
    sFilesOut = sFilesOut{1};
end

end