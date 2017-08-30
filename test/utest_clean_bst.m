function utest_clean_bst()
% Clean brainstorm from anything done during unit test
global GlobalData;
GlobalData.Program.isServer = 0;
GlobalData.Program.HandleExceptionWithBst = 1;

% Delete unit test protocol
ProtocolName = 'nst_utest';
gui_brainstorm('DeleteProtocol', ProtocolName);

% Flush latest error message
GlobalData.lastestFullErrMsg = '';
GlobalData.lastestConsoleMsg = '';
end