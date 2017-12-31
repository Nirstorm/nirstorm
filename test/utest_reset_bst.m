function utest_reset_bst()
% Remove unit test protocol from brainstorm DB
% Reset stored messages.

global GlobalData;

% Delete unit test protocol
ProtocolName = 'nst_utest';
gui_brainstorm('DeleteProtocol', ProtocolName);
db_dir = bst_get('BrainstormDbDir');
nst_protocol_dir = fullfile(db_dir, ProtocolName);
if exist(nst_protocol_dir, 'dir')
    rmdir(nst_protocol_dir, 's');
end

% Flush latest error message
GlobalData.lastestFullErrMsg = '';
GlobalData.lastestConsoleMsg = '';
end
