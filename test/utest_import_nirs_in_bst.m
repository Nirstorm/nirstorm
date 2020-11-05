function sFile = utest_import_nirs_in_bst(nirs_fn, clean, raw_importation)

% RAW=1 -> import as link to raw file
% RAW=0 -> import data into brainstorm DB.

if nargin < 2
    clean = 1; % always clean by default
end

if nargin < 3
    raw_importation = 1;
end

ProtocolName = 'nst_utest';

%% Clean nst_utest protocol if needed
if clean && ~isempty(bst_get('Protocol', ProtocolName))
    % Delete existing protocol
    gui_brainstorm('DeleteProtocol', ProtocolName);
    
    db_dir = bst_get('BrainstormDbDir');
    nst_protocol_dir = fullfile(db_dir, ProtocolName);
    if exist(nst_protocol_dir, 'dir')
        rmdir(nst_protocol_dir, 's');
    end
end

%% Ensure that nst_utest protocol exists
if isempty(bst_get('Protocol', ProtocolName))
    % Create new protocol
    gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
end

%% Resolve file format
[rr, bfn, ext] = fileparts(nirs_fn);
switch(ext)
    case '.nirs' %Homer
        format = 'NIRS-BRS';
    case '.txt' %Artinis TODO: enable importation in Brainstorm
        format = 'NIRS-ARTINIS';   
end

%% Import data as raw or data file
subject_name = 'test_subject';
if raw_importation
    sFile = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                        'subjectname',  subject_name, ...
                        'datafile',       {nirs_fn, format}, ...
                        'channelreplace', 1, ...
                        'channelalign',   0);    
else
    sFile = bst_process('CallProcess', 'process_import_data_time', [], [], ...
                        'subjectname',  subject_name, ...
                        'condition',    'test', ...
                        'datafile',     {nirs_fn, format}, ...
                        'timewindow',   [], ...
                        'split',        0, ...
                        'ignoreshort',  1, ...
                        'channelalign', 0, ...
                        'usectfcomp',   0, ...
                        'usessp',       0, ...
                        'freq',         [], ...
                        'baseline',     []);
end

end
