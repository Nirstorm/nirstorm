function anatomy_based_full_group_pipeline_V1()

 % Example for the template- and surface-based full pipeline, using the
 % function NST_PPL_SURFACE_TEMPLATE_V1.
 %
 % This script downloads some sample data of 10 "dummy" subjects (27 Mb)
 % Total max amount of data to download: (TODO) Mb, the user is asked for download 
 % confirmation.
 %
 % Data are imported in a dedicated protocol, using specific name convention
 % handled by NST_PPL_ANATOMY_BASED_V1.
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
protocol_name = 'TestAnatomyBasedGroupPipelineV1';

if isempty(bst_get('Protocol', protocol_name))
    % Might need a fix
    gui_brainstorm('CreateProtocol', protocol_name, 0, 0);
end
db_save();
 %% Fetch data
subject_names = {'S01',               ...
                 'S04', 'S05',        ...
                 'S07', 'S08', 'S09', ...
                 'S10', 'S11'};
             
nb_subjects = length(subject_names);
% TODO: resolve file names from Subject names
requested_files{4*nb_subjects}={''};

for i=1:nb_subjects
    requested_files{i}={'Tapping',subject_names{i}, 'anatomy' };
    requested_files{i+nb_subjects}={'Tapping',subject_names{i}, 'data', ['subject_' subject_names{i}(2:3) '.nirs'] };
    requested_files{i+2*nb_subjects}={'Tapping',subject_names{i}, 'data', 'optodes.txt' };
    requested_files{i+3*nb_subjects}={'Tapping',subject_names{i}, 'data', 'headpoints' };
end    

data_fns = nst_request_files(requested_files, ...
                             1, nst_get_repository_url(), 1e6, root_folder);
                         
mri_folders=data_fns(1:nb_subjects);
nirs_fns = data_fns((1+nb_subjects):2*nb_subjects);
headpoints=data_fns((1+3*nb_subjects):4*nb_subjects);
%% Import data

options = nst_ppl_surface_template_V1('get_options');

options.import.useDefaultAnat=0;
options.import.mri_folder_type='FreeSurfer';
options.import.nvertices = 2500;
options.import.aseg=1;

options.import.subject(1:nb_subjects)=repmat(options.import.subject,1,nb_subjects);

for i=1:nb_subjects
    options.import.subject{i}.name=subject_names{i};
    options.import.subject{i}.nirs_fn=nirs_fns{i};
    options.import.subject{i}.mri_folder=mri_folders{i};
    options.import.subject{i}.additional_headpoints=headpoints{i};
end    

[sFiles, imported] = nst_ppl_surface_template_V1('import_subjects', options);

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
options.GLM_1st_level.contrasts(1).vector = '[1 0]'; % a vector of weights, as a string 

% Run the pipeline (and  save user markings):
nst_ppl_surface_template_V1('analyse', options, subject_names); % Run the full pipeline
%TODO: full reload of GUI tree at end of pipeline
end