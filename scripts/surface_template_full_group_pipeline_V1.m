function surface_template_full_group_pipeline_V1()

 % Example for the template- and surface-based full pipeline, using the
 % function NST_PPL_SURFACE_TEMPLATE_V1.
 %
 % This script downloads some sample data for 10 "dummy" subjects, as well as
 % the Colin27_4NIRS template if not available (total less than 50 Mb, 
 % the user is asked for download confirmation).
 %
 % Data are imported in a dedicated protocol, using specific name convention
 % handled by NST_PPL_SURFACE_TEMPLATE_V1.
 % Preprocessings and processings are then run until group-level analysis.
 % See NST_PPL_SURFACE_TEMPLATE_V1 for more details.
 %
 % This script illustrates a fully functional analysis pipeline that can
 % serve as a basis for another custom study.
 %
 % For a more detailed description, see the wiki page:
 % TODO: add wiki page
 
 %% Define experiment folder
 % where example data will downloaded and result figures will be stored
 % Default is to create folder in the system temporary directory.
 % Modify root_folder to point to a custom folder.
 % WARNING: it is strongly advised not to use a subdirectory of 
 %          the nirstorm source folder.
 
 tmp_folder = tempdir();
 if ~exist(tmp_folder, 'dir')
     error('Cannot locate temporary folder');
 end
 root_folder = fullfile(tmp_folder, 'nst_ppl_surface_template_V1_example');
 if ~exist(root_folder, 'dir')
    mkdir(root_folder);
 end
 
%% Setup brainstorm
if ~brainstorm('status')
    % Start brainstorm without the GUI if not already running
    brainstorm nogui
end

%% Check Protocol
protocol_name = 'TestSurfaceTemplateGroupPipelineV1';

if isempty(bst_get('Protocol', protocol_name))
    gui_brainstorm('CreateProtocol', protocol_name, 1, 0);
end

% Set template for all subjects
% TODO: make a helper function nst_set_default_template_anatomy()
sTemplates = bst_get('AnatomyDefaults');
iTemplate = strcmpi('Colin27_4NIRS_Jan19', {sTemplates.Name});
if ~any(iTemplate)
    template_bfn = 'Colin27_4NIRS_Jan19.zip';
    template_tmp_fn = nst_request_files({{'template', template_bfn}}, 1, ...
                                        nst_get_repository_url(), 18e6, root_folder);
    % Copy to .brainstorm/defaults/anatomy
    copyfile(template_tmp_fn{1}, ...
             fullfile(bst_get('BrainstormUserDir'), 'defaults', 'anatomy'));
    % Remove temporary file
    rmdir(template_tmp_fn, 's');
end
db_set_template(0, sTemplates(iTemplate), 0);
db_save();
 
 %% Fetch data
subject_names = {'Subject01', 'Subject02', 'Subject03', 'Subject04', ...
                 'Subject05', 'Subject06', 'Subject07', 'Subject08', ...
                 'Subject09', 'Subject10', };
nb_subjects = length(subject_names);
% TODO: resolve file names from Subject names
data_fns = nst_request_files({ {'sample_data', 'template_group_tapping', 'Subject01', 'S01_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject02', 'S02_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject03', 'S03_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject04', 'S04_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject05', 'S05_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject06', 'S06_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject07', 'S07_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject08', 'S08_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject09', 'S09_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject10', 'S10_tapping.nirs'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject01', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject02', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject03', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject04', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject05', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject06', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject07', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject08', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject09', 'optodes.txt'}, ...
                               {'sample_data', 'template_group_tapping', 'Subject10', 'optodes.txt'}, ...
                               }, ...
                             1, nst_get_repository_url(), 1e6, root_folder);
nirs_fns = data_fns(1:nb_subjects);

%% Import data

options = nst_ppl_surface_template_V1('get_options'); % get default pipeline options
 
[sFiles, reimported] = nst_ppl_surface_template_V1('import', options, nirs_fns, ...
                                                   subject_names);
% Read events as extended
for ifile=1:length(sFiles)
    % TODO: use extended events
    % Process: Read from channel
    bst_process('CallProcess', 'process_evt_read', sFiles{ifile}, [], ...
                'stimchan',  'NIRS_AUX', ...
                'trackmode', 3, ...  % Value: detect the changes of channel value
                'zero',      0);
    % Process: Rename event (standard_fix>standard)
    bst_process('CallProcess', 'process_evt_rename', sFiles{ifile}, [], ...
                'src',  'AUX1', ...
                'dest', 'motor');
end

%% Run pipeline
% Can be in another script if manual markings are required

% Specify directories for saving manual markings (not yet implemented)
%   options.moco.export_dir = 'path/to/store/motion_events'
%   options.tag_bad_channels.export_dir = 'path/to/store/bad_channels'

options.GLM_1st_level.contrasts(1).name = 'motor';
options.GLM_1st_level.contrasts(1).vector = [1 0];

% Run the pipeline (and  save user markings):
nst_ppl_surface_template_V1('analyse', options, subject_names); % Run the full pipeline

end