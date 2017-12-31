function sFile = utest_import_nirs_in_bst(nirs_fn)
%% Ensure that nst_utest protocol exists
ProtocolName = 'nst_utest';
% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);

db_dir = bst_get('BrainstormDbDir');
nst_protocol_dir = fullfile(db_dir, ProtocolName);
if exist(nst_protocol_dir, 'dir')
    rmdir(nst_protocol_dir, 's');
end

% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

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
