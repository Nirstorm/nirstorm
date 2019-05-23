function surface_template_full_group_pipeline_V1()
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
 % Data are imported in brainstorm in a dedicated protocol, using specific 
 % naming conventions specific to NST_PPL_SURFACE_TEMPLATE_V1.
 % Preprocessings and processings are then run up to group-level analysis:
 %          1) Resampling to 5hz
 %          2) Bad channels detection
 %          3) Conversion to delta optical density
 %          4) High pass filter
 %          5) Projection on the cortical surface using head model
 %          6) 1st level GLM with pre-coloring
 %          7) group-level GLM with MFX contrast t-maps
 %
 % This script illustrates a minimal fully functional analysis pipeline that can
 % serve as a basis for another custom study.
 %
 % For a more detailed description, see the wiki page:
 % https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-pipeline:-surface-and-template-based
 % 
 % For a more comprehensive example script, see:
 % surface_template_full_group_pipeline_V1_all_opts.m
 
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


%% Import data
options = nst_ppl_surface_template_V1('get_options');

options.import.useDefaultAnat=1;
options.import.subject(1:nb_subjects)=repmat(options.import.subject,1,nb_subjects);

for i=1:nb_subjects
    options.import.subject{i}.name=subject_names{i};
    options.import.subject{i}.nirs_fn=data_fns{i};
    options.import.subject{i}.additional_headpoints=data_fns{i+2*nb_subjects};
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
                    'src',  'AUX1', 'dest', 'motor');
        % Convert to extended event-> add duration of 30 sec to all motor events
        bst_process('CallProcess', 'process_evt_extended', sFiles{ifile}, [], ...
                    'eventname',  'motor', 'timewindow', [0, 30]);
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