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
%   EXTRA can be a string or cell of string corresponding to extra
%   installation scenarios. Files specified in MANIFEST.<extra> will then
%   be installed.
%   For instance 'debug' will install debuging scripts: bst_call.m and bst_process.m.
%         They enable bypassing Brainstorm exception handling.
%         IMPORTANT: 'debug' is needed to run unit tests.
%
%  To cleanly uninstall nirstorm, run nst_uninstall() (see nst_unsintall.m)
%
if nargin < 1           
   mode = 'copy'; 
end

if nargin < 2
    extra = {};
elseif ~iscellstr(extra)
    if ~ischar(extra)
        error('Argument "extra" must be a string');
    end
    extra = {extra};
end


%% Check Brainstorm installation
try
    bst_folder = bst_get('BrainstormUserDir');
catch
    msg = ['Could not find Brainstorm installation. '...
           'Check that matlab path contains Brainstorm folders'];
    throw(MException('Nirstorm:Installation', msg));
end

nistorm_folder      = fileparts(which('nst_install'));
bst_process_folder  = fullfile(bst_folder, 'process');
bst_functions_folder =  bst_fullfile( bst_get('BrainstormUserDir'), 'nirstorm' );


addpath(fullfile(nistorm_folder, 'dist_tools'));
install_package('nirstorm', fullfile(nistorm_folder,'bst_plugin'), bst_process_folder, mode, extra, 0);

% Move the functions to the external folder of brainstorm
file_move(fullfile(bst_get('UserProcessDir'),'nst_*'), bst_functions_folder);
addpath(bst_functions_folder);

% Move the Mex file to the Mex folder
file_move(fullfile(bst_get('UserProcessDir'),'*.mex*'), bst_get('UserMexDir'));


rmpath(fullfile(nistorm_folder, 'dist_tools'));

end