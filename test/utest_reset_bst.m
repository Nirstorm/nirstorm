function utest_reset_bst()
% Remove unit test protocol from brainstorm DB
% Reset stored messages.

global GlobalData;

% Delete unit test protocol
ProtocolName = 'nst_utest';
gui_brainstorm('DeleteProtocol', ProtocolName);

% Flush latest error message
GlobalData.lastestFullErrMsg = '';
GlobalData.lastestConsoleMsg = '';
end