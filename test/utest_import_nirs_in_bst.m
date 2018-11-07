function sFile = utest_import_nirs_in_bst(nirs_fn, clean)

if nargin < 2
    clean = 1; % always clean by default
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

%% Import data as raw file
[rr, bfn, ext] = fileparts(nirs_fn);
switch(ext)
    case '.nirs' %Homer
        sFile = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                            'subjectname',  'test_subject', ...
                            'datafile',       {nirs_fn, 'NIRS-BRS'}, ...
                            'channelreplace', 1, ...
                            'channelalign',   0);
    case '.txt' %Artinis TODO: enable importation in Brainstorm
        sFile = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                            'subjectname',  'test_subject', ...
                            'datafile',       {nirs_fn, 'NIRS-ARTINIS'}, ...
                            'channelreplace', 1, ...
                            'channelalign',   0);   
end
        % [sSubject, iSubject] = bst_get('Subject',  sFile.SubjectName);
end
