function sFilesOut = nst_run_bst_proc(out_items_names, force_redo, ProcessName, sFiles, sFiles2, varargin)
% Allow to run a brainstorm process only if expected outputs do not already exist 
% or if recomputation is forced.
% Outputs have predefined names as given by out_items_names. Each output
% will be unique. This overrides the default behaviour of brainstorm 
% to always run the process and add a suffix to the output name if
% an item with the same name already exists.
% 
% nst_run_bst_proc forces a process to have unique outputs. If outputs
% already exist, the process will not be executed.
% Unless force_redo=1. In this case, the exisiting outputs will be deleted
% prior to executing the process.
% 
% Note that the definition of the output names (Comment field) is
% overriden and no longer determined by the called process. Therefore this
% function is rather suited only for processes with a fixed predictible number 
% of outputs.
%
% Input:
%    - out_items_names (cell array of str):
%          Names of the output items (Comment field). Must be the same length
%          as the actual number of outputs given by process "ProcessName".
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

bst_report('Info', ProcessName, sFiles, 'It is a meeee');
sFilesOut = {}; %stub
end