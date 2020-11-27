function nst_uninstall()
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

nistorm_folder      = fileparts(which('nst_install'));

addpath(fullfile(nistorm_folder, 'dist_tools'));
uninstall_package('nirstorm', bst_process_folder);
rmpath(fullfile(nistorm_folder, 'dist_tools'));

end