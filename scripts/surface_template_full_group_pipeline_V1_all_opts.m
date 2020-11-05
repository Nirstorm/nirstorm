function surface_template_full_group_pipeline_V1_all_opts()
 % Example for the template- and surface-based full pipeline, using the
 % function NST_PPL_SURFACE_TEMPLATE_V1.
 %
 % This script downloads some sample data of 10 "dummy" subjects (27 Mb), 
 % as well as the Colin27_4NIRS template (19 Mb) if not available.
 % For the analysis part, precomputed fluence data are also downloaded.
 % Total max amount of data to download: 50 Mb, the user is asked for download 
 % confirmation.
 % All downloaded data will be stored in .brainstorm/defaults/nirstorm
 %
 % Since several files will be written outside of the brainstorm database,
 % a temporary folder will be used (see output_dir).
 %
 % Data are imported in brainstorm in a dedicated protocol, using specific 
 % naming conventions handled by NST_PPL_SURFACE_TEMPLATE_V1.
 % Then, this script emulates user-defined movement artefacts markings as
 % well as bad channel taggings.
 % These markings will be saved in output_dir. Everytime NST_PPL_SURFACE_TEMPLATE_V1
 % is run, markings in the brainstorm DB are saved in the output_dir.
 % When some nirs data is not available in the DB and need (re)importation,
 % these markings are loaded at the same time.
 %
 % Preprocessings and processings are then run up to group-level analysis:
 %          0) Figure outputs of raw time-series and SCI maps
 %          1) Motion-correction
 %          2) Resampling to 5hz
 %          3) Bad channels detection
 %          4) Conversion to delta optical density
 %          5) High pass filter
 %          6) Figure output of preprocessed time-series
 %          5) Projection on the cortical surface using head model
 %          6) 1st level GLM with pre-coloring
 %          7) Subject-level t-maps with figure outputs
 %          8) group-level GLM with MFX contrast t-maps
 %
 % This script illustrates a minimal fully functional analysis pipeline that can
 % serve as a basis for another custom study.
 %
 % For a more detailed description, see the wiki page:
 % https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-pipeline:-surface-and-template-based#example-script-with-all-options
 % 
 % For a minimal example script, see:
 % surface_template_full_group_pipeline_V1.m
 %
 
%% Setup brainstorm
if ~brainstorm('status')
    % Start brainstorm without the GUI if not already running
    brainstorm nogui
end

%% Check Protocol
protocol_name = 'TestSurfaceTemplateGroupPipelineV1';
if isempty(bst_get('Protocol', protocol_name))
    gui_brainstorm('CreateProtocol', protocol_name, 1, 0); % UseDefaultAnat=1, UseDefaultChannel=0
end

% Set template for default anatomy
nst_bst_set_template_anatomy('Colin27_4NIRS_Jan19');
 
%% Fetch data
% Get list of local nirs files for the group data.
% The function nst_io_fetch_sample_data takes care of downloading data to
% .brainstorm/defaults/nirstorm/sample_data if necessary
[nirs_fns, subject_names] = nst_io_fetch_sample_data('template_group_tapping'); 

%% Set root output dir (to store markings, figures and results)
tmp_folder = tempdir();
if ~exist(tmp_folder, 'dir')
    error('Cannot locate temporary folder');
end
output_dir = fullfile(tmp_folder, 'nst_ppl_surface_template_V1_example');
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end

%% Get default options
options = nst_ppl_surface_template_V1('get_options'); % get default pipeline options 

%% Set output directories for markings
options.moco.export_dir = fullfile(output_dir, 'markings_moco');
options.tag_bad_channels.export_dir = fullfile(output_dir, 'markings_bad_channels');

%% Import data
sFiles = nst_ppl_surface_template_V1('import', options, nirs_fns, subject_names);

% Read stimulation events from AUX channel
for ifile=1:length(sFiles)
    evt_data = load(file_fullpath(sFiles{ifile}), 'Events');
    if ~any(strcmp({evt_data.Events.label}, 'motor')) % Insure that events were not already loaded
        bst_process('CallProcess', 'process_evt_read', sFiles{ifile}, [], ...
                    'stimchan',  'NIRS_AUX', ...
                    'trackmode', 3, ...  % Value: detect the changes of channel value
                    'zero',      0);
        % Rename event AUX1 -> motor
        bst_process('CallProcess', 'process_evt_rename', sFiles{ifile}, [], ...
                    'src',  'AUX1', 'dest', 'motor');
        % Convert to extended event-> add duration of 30 sec to all motor events
        bst_process('CallProcess', 'process_evt_extended', sFiles{ifile}, [], ...
                    'eventname',  'motor', 'timewindow', [0, 30]);
    end
end

%% Simulate user-inputs
% Add tagging of movement artefacts

% Add bad channel tagging

%% Run pipeline
options.GLM_1st_level.stimulation_events = {'motor'};
options.GLM_1st_level.contrasts(1).label = 'motor';
options.GLM_1st_level.contrasts(1).vector = '[1 0]'; % a vector of weights, as a string 

% Run the pipeline (and  save user markings):
nst_ppl_surface_template_V1('analyse', options, subject_names); % Run the full pipeline
%TODO: full reload of GUI tree at end of pipeline
end
