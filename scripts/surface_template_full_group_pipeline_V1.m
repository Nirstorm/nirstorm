function surface_template_full_group_pipeline_V1()

 % Example for the template- and surface-based full pipeline, using the
 % function NST_PPL_SURFACE_TEMPLATE_V1.
 %
 % This script downloads some sample data of 10 "dummy" subjects (27 Mb), 
 % as well as the Colin27_4NIRS template (19 Mb) if not available.
 % For the analysis part, precomputed fluence data are also downloaded.
 % Total max amount of data to download: 50 Mb, the user is asked for download 
 % confirmation.
 %
 % Data are imported in a dedicated protocol, using specific name convention
 % handled by NST_PPL_SURFACE_TEMPLATE_V1.
 % Preprocessings and processings are then run until group-level analysis:
 %          1) Resampling to 5hz
 %          2) Detect bad channels
 %          3) Convert to delta optical density
 %          4) High pass filter
 %          5) Project on the cortical surface using head
 %          6) 1st level GLM (pre-coloring)
 %          7 level GLM with MFX contrast t-maps
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
[sFiles, imported] = nst_ppl_surface_template_V1('import', options, nirs_fns, subject_names);

% Read stimulation events from AUX channel
for ifile=1:length(sFiles)
    if imported(ifile)
        % Read events from aux channel
        bst_process('CallProcess', 'process_evt_read', sFiles{ifile}, [], ...
                    'stimchan',  'NIRS_AUX', ...
                    'trackmode', 3, ...  % Value: detect the changes of channel value
                    'zero',      0);
        % Rename event AUX1 -> motor
        bst_process('CallProcess', 'process_evt_rename', sFiles{ifile}, [], ...
                    'src',  'AUX1', ...
                    'dest', 'motor');
        % Convert to extended event-> add duration of 30 sec to all motor events
        bst_process('CallProcess', 'process_evt_extended', sFiles{ifile}, [], ...
                    'eventname',  'motor', ...
                    'timewindow', [0, 30]);
    end
end

%% Run pipeline
options.GLM_1st_level.stimulation_events = {'motor'};
options.GLM_1st_level.contrasts(1).label = 'motor';
options.GLM_1st_level.contrasts(1).vector = '[1 0]'; % a string

% Run the pipeline (and  save user markings):
nst_ppl_surface_template_V1('analyse', options, subject_names); % Run the full pipeline
%TODO: full reload of GUI tree at end of pipeline
end