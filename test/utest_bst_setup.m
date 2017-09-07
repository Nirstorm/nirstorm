function utest_bst_setup()
% Run brainstorm with no user interface, if not already running.
% Force brainstorm to run in server mode (no user interaction).
% Activate exception bypassing. It requires to override bst_error, bst_call
% and bst_process with the ones from nirstorm 
% -> use nst_install('link', 'debug') or nst_install('copy', 'debug')
%
global GlobalData;
if ~brainstorm('status')
    brainstorm nogui;
end
GlobalData.Program.isServer = 1;
GlobalData.Program.HandleExceptionWithBst = 0;
GlobalData.lastestFullErrMsg = '';
GlobalData.lastestConsoleMsg = '';
panel_process_select('ParseProcessFolder', 1);
end
