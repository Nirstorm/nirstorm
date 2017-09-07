function nst_install(mode, extra)
%NST_INSTALL Installation of NIRSTORM plugin for Brainstorm (linux and windows).
%   NST_INSTALL(MODE) install NIRSTORM processes to brainstorm user folder:
%     - linux: $HOME/.brainstorm/process
%     - windows: C:\Documents and Settings\username\.brainstorm
%   It backups existing scripts in brainstorm user folder.
%
%   IMPORTANT: assume brainstorm has been installed and its functions are
%   available in matlab's path.
%
%   If MODE == 'copy' then all processes are copied to the installation
%   folder. This is recommended for an end-user installation (no frequent
%   code updates). 
%   If MODE == 'link' (linux only) then symbolink links pointing to scripts 
%   in the  source folder are created in the installation folder. 
%   This is recommeded when source code needs to be updated often. 
%   Any modification in the source folder will be available in brainstorm.
%   IMPORTANT: for *new scripts*, installation has to be run again to create
%   new symbolic links.

%   NST_INSTALL(MODE, EXTRA)
%   Install extra scripts that override some brainstorm functions.
%   WARNING: these are mostly in-dev features so they may be unstable.
%            They are strongly dependent on the current version of
%            brainstorm, so it's better to have the most up-to-date
%            version.
%            DO NOT install these scripts unless you know what you're doing ;)
%   EXTRA can be one of:
%       - 'none' (default): do not install extra scripts
%       - 'debug': install debuging scripts: bst_call.m and bst_process.m.
%            They enable bypassing Brainstorm exception handling.
%         IMPORTANT: this is needed to run unit tests.
%
%  To cleanly uninstall nirstorm, run nst_uninstall() (see nst_unsintall.m)
%
if nargin < 1           
   mode = 'copy'; 
end

if nargin < 2
    extra = {};
else
    extra = {extra};
end

if nargin < 4
    dry = 0;
end

%% Check Brainstorm installation
try
    bst_folder = bst_get('BrainstormUserDir');
catch
    msg = ['Could not find Brainstorm installation. '...
           'Check that matlab path contains Brainstorm folders'];
    throw(MException('Nirstorm:Installation', msg));
end

bst_process_folder = fullfile(bst_folder, 'process');
if ~exist(bst_process_folder, 'dir')
    display(['Could not find Brainstorm process folder "' ...
             bst_process_folder '". Check brainstorm installation']);
    return;
end
addpath(fullfile(pwd, 'dist_tools'));
install_package('nirstorm', 'bst_plugin', bst_process_folder, mode, extra, dry);
