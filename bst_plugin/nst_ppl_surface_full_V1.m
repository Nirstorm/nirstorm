function varargout = nst_ppl_surface_full_V1(action, options, arg1, arg2)
%NST_PPL_SURFACE_TEMPLATE_V1
% Manage a full template- and surface-based pipeline starting from raw NIRS data
% up to GLM group analysis (if enough subjects).
%
% IMPORTANT: although each subject has its own optode coordinate file during importation,
% the same optode coordinates are used for all subjects. These coordinates
% are the one from the first given subject.
%
% This pipeline can keep track of user-defined markings outside of brainstorm db 
% such as movement events and bad channels. This allows to safely flush all
% brainstorm data while keeping markings.
% 
% This function is intended to be called from batch scripts where the user
% can add some custom steps. Here is the workflow:
%
%   options = NST_PPL_SURFACE_TEMPLATE_V1('get_options'); % get default pipeline options
%
%   % Define import options (optional):
%   options.moco.export_dir = 'path/to/store/motion_events'
%   options.tag_bad_channels.export_dir = 'path/to/store/bad_channels'
%
%   % Import some nirs data along with event markings:
%   subject_names = {'subj1', 'subj2'};
%   sFilesRaw = NST_PPL_SURFACE_TEMPLATE_V1('import', options, {'data1.nirs', 'data2.nirs'}, subject_names);
%   for ifile=1:length(sFilesRaw)
%     % Tweak sFilesRaw{ifile} here, eg import stimulation event.
%   end
%
%   % User can manually tag motion events and bad channels here
%
%   % Customize options:
%   options.GLM_1st_level.contrasts(1).name = 'my_contrast1';
%   options.GLM_1st_level.contrasts(1).vector = [0 1 -1 0];
% 
%   % Run the pipeline (and  save user markings):
%   NST_PPL_SURFACE_TEMPLATE_V1('analyse', options, subject_names); % Run the full pipeline
%
%   % For a working example see
%   nirstorm/script/surface_template_full_group_pipeline.m
%
% DEFAULT_OPTIONS = NST_PPL_SURFACE_TEMPLATE_V1('get_options')
%     Return default options
%
% FILES_RAW = NST_PPL_SURFACE_TEMPLATE_V1('import', OPTIONS, NIRS_FNS, SUBJECT_NAMES)
%     Import all nirs files in database and use given subjects (skip if exists).
%     NIRS_FNS is a cell array of str.
%     If SUBJECT_NAMES is empty or not given, then use base filename as
%     subject names. If not empty, then it must be a cell array of str with the 
%     same length as NIRS_FNS.
%
%     Used options:
%        - options.import.redo
%
%     Return:
%         FILES_RAW: brainstorm file pathes to imported data.
%
%  NST_PPL_SURFACE_TEMPLATE_V1('analyse', OPTIONS, GROUPS | SUBJECT_NAMES)
%   
%     Apply pipeline to given group(s) of subjects.
%     ASSUME: all subjects in a given protocol have the same template
%             anatomy. See NST_PPL_SURFACE_TEMPLATE_V1('import')
%
%     List of steps:
%        - For 1st subject: head model for all pairs
%            Clone first item of FILES_RAW and precompute a head model for
%            all possible pairs, which can be used for other subjects.
%            -> head models are not recomputed for each subject.
%            Create a dummy subject called "full_head_model...".
%        Then for each subject:
%          -) Export user-defined inputs [optional]: TODO
%             - movement events
%             - bad channels
%          1) Motion correction
%              ASSUME: event group "movement_artefacts" exists in each FILES_RAW 
%                      and has been filled by user before calling this function.
%                      Note that NST_PPL_SURFACE_TEMPLATE_V1('import')
%                      creates this event group if necessary.
%          2) Resampling:
%              TODO: check interpolation errors when there are spikes
%          3) Detect bad channels
%          4) Convert to delta optical density
%          5) High pass filter
%          6) Compute head model (from "full_head_model...")
%          7) Project on the cortical surface
%          8) 1st level GLM:
%              - build design matrix from stimulation events
%              - OLS fit with pre-coloring
%              - compute contrasts
%        Group-level:
%          1) 2nd level GLM:
%              - build design matrix  
%              - OLS fit
%              - MFX contrast t-maps
%          2) Extract group-masked subject-level maps [optional]
%
% TODO:
% - handle when no contrast defined
% - options documentation
% - export manual inputs
% - importation of manual inputs:
% - wiki page
% - utest
% 
global GlobalData;

assert(ischar(action));

%TODO check options when given

switch action
    case 'get_options'
        if nargin > 1
            protocol_name = options;
            assert(ischar(protocol_name));
            varargout{1} = get_options(protocol_name);
        else
            varargout{1} = get_options();
        end
        return;
        
    case 'rois_summary.get_mask_combinations'
        varargout{1} = get_mask_combinations();
        return;
    case 'import_nirs'
        if nargin >= 4
            subject_names = arg2;
        else
            subject_names = cell(size(arg1));
            subject_names(:) = {''};
        end
        [imported_files, redone] = import_nirs_files(arg1, subject_names, options);
        varargout{1} = imported_files;
        varargout{2} = redone;
        return;
    case 'import_subjects'
        [imported_files, redone] = import_subjects(options);
        
        varargout{1} = imported_files;
        varargout{2} = redone;
        return;
    case 'preprocessing'
        if nargin >= 3
            % Only compute the fluences for the specified subjects
            subject_names = arg1;
        else
            % Compute the fluences for every subjects 
            subject_names = cell(size(options.import.subject));
            for iSubject=1:size(options.import.subject,2)
                subject_names{iSubject} = options.import.subject{iSubject}.name;
            end    
        end
         sFiles=cell(size(options.import.subject));
         redones=cell(size(options.import.subject));
         
        for iSubject=1:length(subject_names)
            [sFiles{iSubject},hb_types,redones{iSubject}, preproc_folder] = preprocs(subject_names{iSubject},options);
            varargout{1} = sFiles;
            varargout{2} = redones;
        end
        return;
    case 'compute_fluences'
        if nargin >= 3
            % Only compute the fluences for the specified subjects
            subject_names = arg1;
        else
            % Compute the fluences for every subjects 
            subject_names = cell(size(options.import.subject));
            for iSubject=1:size(options.import.subject,2)
                subject_names{iSubject} = options.import.subject{iSubject}.name;
            end    
        end
       
        compute_fluences(subject_names,options);
        return
    case 'analyse'
    otherwise
        error('Unknown action: %s', action);
end

if isempty(arg1)
    error('Empty input group or subjects definition');
else
    if isstruct(arg1)
        groups = arg1;
        assert(isfield(groups, 'label'));
        assert(isfield(groups, 'subject_names'));
    else
        groups.label = '';
        groups.subject_names = arg1;
    end
    %TODO: check groups -> warn if they are overlapping
end

if strcmp(options.save_fig_method, 'export_fig') && ~function_exists('export_fig')
    error('"export_fig" not found. Can be installed from "https://github.com/altmany/export_fig"');
end

force_redo = options.redo_all;

create_dir(options.fig_dir);
create_dir(options.moco.export_dir);
create_dir(options.tag_bad_channels.export_dir);
create_dir(options.GLM_group.rois_summary.csv_export_output_dir);


nb_groups = length(groups);
group_condition_names = cell(1, nb_groups);
any_rois_summary_redone = 0;
all_sFiles_con = cell(1, nb_groups);
for igroup=1:nb_groups
    redo_group = options.GLM_group.redo;
    subject_names = groups(igroup).subject_names;
    group_label = groups(igroup).label;
    
    %% Within-subject analyses
    for isubject=1:length(subject_names)

        subject_name = subject_names{isubject};


        
        % Run preprocessings
        [sFiles_preprocessed, hb_types, redone_preprocs, preproc_folder] = preprocs(subject_name,options, force_redo);
        
        % Run 1st level GLM
        [sFiles_GLM, sFiles_con, redone_any_contrast, glm_folder] = glm_1st_level(sFiles_preprocessed, options, ...
                                                                                  redone_preprocs | force_redo);
        
        if isubject==1
            nb_hb_types = size(sFiles_con, 1);
            nb_contrasts = size(sFiles_con, 2);
            all_sFiles_con{igroup} = cell(nb_hb_types, nb_contrasts, length(subject_names));
        else
            assert(size(sFiles_con, 1) == nb_hb_types);
            assert(size(sFiles_con, 2) == nb_contrasts);
        end
        all_sFiles_con{igroup}(:,:, isubject) = sFiles_con;

        if options.clean_preprocessings
            full_preproc_folder = fileparts(sFiles_preprocessed{1});
            [sStudy, iStudy] = bst_get('StudyWithCondition', full_preproc_folder);
            db_delete_studies(iStudy);
        end

        redo_group = redo_group | redone_any_contrast;
    end



%% Finalize
if prev_iCortex ~= iCortex
   % Set default cortical surface to original one
   db_surface_default(0, 'Cortex', prev_iCortex);
   panel_protocols('RepaintTree');
end
end
end

function [sFilesHbProj, hb_types, redone, preproc_folder]  = preprocs(subject_name,options, force_redo)

if nargin < 4
    force_redo = 0;
end

sFile_raw = nst_get_bst_func_files(subject_name, ['origin' get_ppl_tag()], 'Raw');
if isempty(sFile_raw)
    error(sprintf('Cannot find "origin/Raw" data for subject "%s". Consider using nst_ppl_surface_template_V1(''import'',...).',  ...
                  subject_name)); %#ok<SPERR>
end
preproc_folder = sprintf('preprocessing%s/', get_ppl_tag());

% Compute Scalp coupling index
nst_run_bst_proc([preproc_folder 'SCI'], force_redo | options.sci.redo, 'process_nst_sci', sFile_raw);

% TODO: export motion correction tagging to external file
% sRaw = load(file_fullpath(sFile_raw));
% sExport.Events = sRaw.Events(strcmp({sRaw.Events.label}, 'movement_artefacts'));
% export_events(sExport, [], moco_export_fn);

% TODO: export bad channel tagging information

% TODO: plot raw input signals
% fig_bfn = sprintf('%s_%s_signals_raw.png', SubjectName, data_tag);
% fig_fn = protect_fn_str(fullfile(options.fig_dir, fig_bfn ));
% if ~isempty(options.fig_dir) && options.make_figs && options.plot_raw_signals.do && ...
%         (force_redo || options.plot_raw_signals.redo || ~exist(fig_fn, 'file'))
%    plot_signals(sFile_raw, fig_fn, options);
% end

% Deglitching
if options.deglitch.do
    redo_parent = force_redo | options.deglitch.redo;
    sFile_deglitched = nst_run_bst_proc([preproc_folder 'Deglitched'], redo_parent, ...
                                        'process_nst_deglitch', sFile_raw, [], ...
                                        'factor_std_grad', options.deglitch.agrad_std_factor);
else
    redo_parent = force_redo;
    sFile_deglitched = sFile_raw;
end

% Motion correction
redo_parent = redo_parent | options.moco.redo;
[sFileMoco, redo_parent] = nst_run_bst_proc([preproc_folder 'Motion-corrected'], redo_parent, ...
                             'process_nst_motion_correction', sFile_deglitched, [], ...
                             'option_event_name', 'movement_artefacts');
                         
% Resample to 5Hz (save some space)
redo_parent = redo_parent | options.resample.redo;
[sFileMocoResampled, redo_parent] = nst_run_bst_proc([preproc_folder 'Motion-corrected | Resampled'], redo_parent, ...
                                                      'process_resample', sFileMoco, [], ...
                                                      'freq', options.resample.freq, ...
                                                      'read_all', 1);
% Process: Detect bad channels
% This one is done in-place -> not tracked to handle do/redo scenarios
if redo_parent
    bst_process('CallProcess', 'process_nst_detect_bad', sFileMocoResampled, [], ...
                'option_remove_negative', 1, ...
                'option_invalidate_paired_channels', 1, ...
                'option_max_sat_prop', options.tag_bad_channels.max_prop_sat_ceil, ...
                'option_min_sat_prop', options.tag_bad_channels.max_prop_sat_floor);
end
% Convert to delta OD
redo_parent = redo_parent | options.dOD.redo;
[sFile_dOD, redo_parent] = nst_run_bst_proc([preproc_folder 'dOD'], redo_parent, 'process_nst_dOD', sFileMocoResampled, [], ...
                                            'option_baseline_method', options.dOD.baseline_def); 
              
% Band pass filter
redo_parent = redo_parent | options.high_pass_filter.redo;
[sFile_dOD_filtered, redo_parent] = nst_run_bst_proc([preproc_folder 'dOD | filtered'],  redo_parent, 'process_bandpass', sFile_dOD, [], ...
                                      'highpass', {options.high_pass_filter.low_cutoff, ''}, ...
                                      'lowpass', {0, ''}, ...
                                      'attenuation', 'relax', ...
                                      'mirror', 0, ...
                                      'sensortypes', 'NIRS');
                                  
                                  
% Compute head model from full head model
redo_parent = redo_parent | options.head_model.redo;
fluences_dir=fullfile(options.fluences.export_dir,subject_name);

[dummy_out, redo_parent] = nst_run_bst_proc([preproc_folder 'head model'], redo_parent, ...
                           'process_nst_import_head_model', sFile_dOD_filtered, [], ...
                           'data_source', fluences_dir ,...
                           'use_closest_wl', 1, 'use_all_pairs', 0, ...
                           'force_median_spread', 0, ...
                           'normalize_fluence', 1, ...
                           'smoothing_fwhm', options.head_model.smoothing_fwhm);
    
                                      


% Project and convert to d[HbX]
redo_parent = redo_parent | options.projection.redo;
proj_method =  options.projection.method;
[sFilesHbProj, redo_parent] = nst_run_bst_proc({[preproc_folder 'dHbO_cortex'], [preproc_folder 'dHbR_cortex']},  redo_parent, ... 
                                                'process_nst_cortical_projection', sFile_dOD_filtered, [], ...
                                                'method', proj_method, ...
                                                'sparse_storage', options.projection.sparse_storage);
hb_types = process_nst_cortical_projection('get_hb_types');
redone = redo_parent;                                  
end

function compute_fluences(subject_names,options)
% Todo : make sure that each subject has been proprely imported 
% Check if the fluences have already been computed and re-computed them
% only if asked (flag redo)

create_dir(options.fluences.export_dir);
outputdir={};
for iSubject=1:length(subject_names)
    
    file_raw = nst_get_bst_func_files(subject_names{iSubject}, ['origin' get_ppl_tag()], 'Raw');
    if isempty(file_raw)
        error(sprintf('Cannot find "origin/Raw" data for subject "%s". Consider using nst_ppl_surface_template_V1(''import'',...).',  ...
                      subject_names{iSubject})); %#ok<SPERR>
    end
        
    outputdir{1}=fullfile(options.fluences.export_dir,subject_names{iSubject});
    create_dir(outputdir{1})
    
    bst_process('CallProcess', 'process_nst_ComputeFluencesforOptodes', file_raw, [], ...
                'outputdir', outputdir, ...
                'mcxlab_nphoton', options.fluences.nphoton, ...
                'mcxlab_flag_thresh',options.fluences.thresh, ...
                'vox2ras',     1);

end
end

function [sFiles_GLM, sFiles_con, redone_any_contrast, glm_folder] = glm_1st_level(sFiles, options, force_redo)

if nargin < 3
    force_redo = 0;
end

stim_events = options.GLM_1st_level.stimulation_events;
if isempty(stim_events)
   error('Stimulation events not defined for building the design matrix.');
   %TODO: use all found events?
end

contrasts = options.GLM_1st_level.contrasts;
if isempty(contrasts)
    warning('Contrasts not defined. Using default single-condition contrasts.');
    contrasts = nst_make_basic_contrasts(stim_events);
end

glm_folder = sprintf('GLM%s/', get_ppl_tag());

[SubjectName, preprocs_folder] = bst_fileparts(bst_fileparts(sFiles{1}), 1);
sSubject = bst_get('Subject', SubjectName);

redone_any_contrast = 0; % Track if any contrast for any file had to be recomputed
sFiles_GLM = cell(1, length(sFiles));

redo_parent = force_redo | options.GLM_1st_level.redo;
for ifile=1:length(sFiles)
    data_cmt = load(file_fullpath(sFiles{ifile}), 'Comment');
    
    comment_glm_prefix{ifile} = [glm_folder 'GLM ' data_cmt.Comment];
    % Process: GLM - design and fit
    [sFiles_GLM{ifile}, redone_fit] = nst_run_bst_proc([comment_glm_prefix{ifile} ' | fitted model'], redo_parent, ...
                                                       'process_nst_glm_fit', sFiles{ifile}, [], ...
                                                       'stim_events',    strjoin(stim_events, ', '), ...
                                                       'hrf_model',      1, ...  % CANONICAL
                                                       'trend',          1, ...
                                                       'fitting',        1, ...  % OLS - precoloring
                                                       'hpf_low_cutoff', options.high_pass_filter.low_cutoff, ...
                                                       'trim_start', options.GLM_1st_level.trim_start, ...
                                                       'save_residuals', 0, ...
                                                       'save_betas',     0, ...
                                                       'save_fit',       0);
end

sFiles_con = cell(length(sFiles_GLM), length(contrasts));
for ifile=1:length(sFiles_GLM)
    for icon=1:length(contrasts)
        % Process: GLM - intra subject contrast
        redo = redone_fit | options.GLM_1st_level.contrast_redo;
        
        [sFile_GLM_con, redone_con] = nst_run_bst_proc([comment_glm_prefix{ifile} ' | con ' contrasts(icon).label], redo, ...
                                                       'process_nst_glm_contrast', sFiles_GLM{ifile}, [], ...
                                                       'Contrast', contrasts(icon).vector);
        sFiles_con{ifile,icon} = sFile_GLM_con;
        
        % GLM - tmaps
        if options.GLM_1st_level.contrast_tstat.do
            redo = redone_con | options.GLM_1st_level.contrast_tstat.redo;
            sFile_GLM_ttest = nst_run_bst_proc([comment_glm_prefix{ifile} ' | con_t-+ ' contrasts(icon).label], redo, ...
                                               'process_nst_glm_contrast_ttest', sFile_GLM_con, [], ...
                                               'tail', 'two');
        
            % Plots
            data_tag = get_bst_file_tag(sFiles{ifile});
            fig_bfn = sprintf('%s_%s_tmap_mcc_%s_pv_thresh_%s_%s.png', ...
                SubjectName, data_tag, options.GLM_1st_level.contrast_tstat.plot.pvalue_mcc_method,...
                nst_format_pval(options.GLM_1st_level.contrast_tstat.plot.pvalue_threshold), ...
                contrasts(icon).label);
            fig_fn = protect_fn_str(fullfile(options.fig_dir, fig_bfn ));
            if ~isempty(options.fig_dir) && options.make_figs && ...
                    options.GLM_1st_level.contrast_tstat.plot.do && ...
                    (redo || options.GLM_1st_level.contrast_tstat.plot.redo || ~exist(fig_fn, 'file'))
                hFigSurfData = view_surface_data(sSubject.Surface(sSubject.iCortex).FileName, ...
                    sFile_GLM_ttest, 'NIRS', 'NewFigure');
                StatThreshOptions = bst_get('StatThreshOptions');
                StatThreshOptions.pThreshold = options.GLM_1st_level.contrast_tstat.plot.pvalue_threshold;
                StatThreshOptions.Correction = options.GLM_1st_level.contrast_tstat.plot.pvalue_mcc_method;
                %StatThreshOptions.Control    = [1 2 3]; % ???
                bst_set('StatThreshOptions', StatThreshOptions);
            
                % TODO: set surface smoothing
                % TODO: set better colormap that does not span values betwn -3 and 3
                if ~isempty(options.fig_background)
                    bst_figures('SetBackgroundColor', hFigSurfData, options.fig_background);
                end
                bst_colormaps('SetDisplayColorbar', 'stat2', 0);
                view(options.fig_cortex_view);
                zoom(hFigSurfData, options.fig_cortex_zoom);
                nst_save_figure(fig_fn, options, hFigSurfData);
                close(hFigSurfData);
            end
        end
        
        redone_any_contrast = redone_any_contrast | redone_con;
    end
end
end


function options = get_options()

options.redo_all = 0;


options.import.redo = 0;
options.import.useDefaultAnat = 1;
options.import.mri_folder_type='FreeSurfer';
options.import.nvertices = 2500;
options.import.aseg=1;

options.import.subject{1}.name='';
options.import.subject{1}.nirs_fn='';
options.import.subject{1}.mri_folder='';
options.import.subject{1}.additional_headpoints='';



options.head_model.surface = 'cortex_lowres';
options.sci.redo = 0;
options.head_model.redo = 0;
options.head_model.smoothing_fwhm=0;


options.fluences.export_dir=''; % where to export the fluences
                                % will be used to compute the head model
                                % 
options.fluences.nphoton=100; % Numner of photon used for the simulation in Million
options.fluences.thresh=0; % between 0 and 1  

options.deglitch.do = 0;
options.deglitch.redo = 0;
options.deglitch.agrad_std_factor = 2.5;

options.moco.redo = 0;
options.moco.export_dir = ''; % Where to export motion correction events manually tagged (for backup) 
                              % -> will be exported before running the analysis
                              % -> will be reimported everytime the importation
                              %    stage is run
options.resample.redo = 0;
options.resample.freq = 5; % Hz

options.dOD.redo = 0;
options.dOD.baseline_def = 1; % 1: mean, 2: median

options.high_pass_filter.redo = 0;
options.high_pass_filter.low_cutoff = 0.01; %Hz

options.tag_bad_channels.redo = 0;
options.tag_bad_channels.max_prop_sat_ceil = 1; % no tagging
options.tag_bad_channels.max_prop_sat_floor = 1; % no tagging
options.tag_bad_channels.export_dir = ''; % Where to export bad channel taggings (backup) 
                                          % -> will be exported before
                                          %    running the analysis
                                          % -> will be reimported everytime 
                                          %    the importation stage is run 

options.projection.redo = 0;
proj_methods = process_nst_cortical_projection('methods');
options.projection.method = proj_methods.Sensitivity_based_interpolation; % proj_methods.MNE;

options.projection.sparse_storage = 0;

options.clean_preprocessings = 0;

options.GLM_1st_level.redo = 0;
options.GLM_1st_level.trim_start = 0; % sec

options.GLM_1st_level.contrast_redo = 0;
options.GLM_1st_level.stimulation_events = [];
options.GLM_1st_level.contrasts = []; 

options.GLM_1st_level.contrast_tstat.do = 0; % not active by default -> only beta values are mandatory for group-level analysis
options.GLM_1st_level.contrast_tstat.redo = 0;
options.GLM_1st_level.contrast_tstat.plot.do = 0; % not active by default -> can produce a lot of figures
options.GLM_1st_level.contrast_tstat.plot.redo = 0; 
options.GLM_1st_level.contrast_tstat.plot.pvalue_threshold = 0.05;
options.GLM_1st_level.contrast_tstat.plot.pvalue_mcc_method = 'none';

options.GLM_group.do = 1;
options.GLM_group.redo = 0;

options.GLM_group.contrast_tstat.plot.do = 1;
options.GLM_group.contrast_tstat.plot.redo = 0;
options.GLM_group.contrast_tstat.plot.pvalue_threshold = 0.001;
options.GLM_group.contrast_tstat.plot.pvalue_mcc_method = 'none';

options.GLM_group.rois_summary.do = 0;
options.GLM_group.rois_summary.atlas = 'MarsAtlas';
options.GLM_group.rois_summary.matrix_col_prefix = '';
options.GLM_group.rois_summary.csv_export_output_dir = '';
mask_combinations = get_mask_combinations();
options.GLM_group.rois_summary.group_masks_combination = mask_combinations.none; % mask_combinations.intersection, mask_combinations.union
options.GLM_group.rois_summary.group_masks_combine_contrasts = 0; 

options.make_figs = 1;
options.save_fig_method = 'saveas'; % 'saveas', 'export_fig'
options.export_fig_dpi = 90;
options.fig_dir = '';
options.fig_background = []; % use default background
options.fig_cortex_view = [89 -24]; % Azimuth and Elevation
                                    % to adjust them manually, right-click on fig 
                                    % then Figure > Matlab controls
                                    % use rotate 3D tool, while moving Az and El
                                    % are displayed in the bottom right of the figure.
options.fig_cortex_smooth = 0;
options.fig_cortex_show_sulci = 0;
options.fig_cortex_zoom = 1;

end

function ptag = get_ppl_tag()
ptag = '__nspst_V1';
end

function operations = get_mask_combinations()
operations = process_nst_combine_masks('get_mask_combinations');
operations.none = length(fieldnames(operations)) + 1;
end


function [files_in, redone_imports] = import_nirs_files(nirs_fns, subject_names, options)
files_in = cell(size(nirs_fns));
redone_imports = zeros(size(nirs_fns));
for ifile=1:length(nirs_fns)
    %% Import data
    nirs_fn = nirs_fns{ifile};
    if isempty(subject_names{ifile})
        [root, subject_name, ext] = fileparts(nirs_fn);
    else
        subject_name = subject_names{ifile};
    end
    condition = ['origin' get_ppl_tag()];
    [file_in, redone] = nst_run_bst_proc([subject_name '/' condition '/Raw'], options.import.redo, ...
                                           'process_import_data_time', [], [], ...
                                           'subjectname',  subject_name, ...
                                           'condition',    condition, ...
                                           'datafile',     {nirs_fn, 'NIRS-BRS'}, ...
                                           'timewindow',   [], ...
                                           'split',        0, ...
                                           'ignoreshort',  1, ...
                                           'channelalign', 1, ...
                                           'usectfcomp',   0, ...
                                           'usessp',       0, ...
                                           'freq',         [], ...
                                           'baseline',     []);
    redone_imports(ifile) = redone;
    %% Manage movement event markings TODO
    if redone
        evt_formats = bst_get('FileFilters', 'events');
        evt_format = evt_formats(strcmp('BST', evt_formats(:,3)), :);

        moco_fn = get_moco_markings_fn(subject_name, options.moco.export_dir);
        if exist(moco_fn, 'file')
            % Load event from pre-saved file
            % TODO: test
            sFile_in = load(file_fullpath(file_in));
            [sFile_in, events] = import_events(sFile_in, [], moco_fn, evt_format);   
        else
            % Create empty event group
            movement_events = db_template('event');
            movement_events.label = 'movement_artefacts';
            sFile_in = bst_process('GetInputStruct', file_in);
            process_nst_import_csv_events('import_events', [], sFile_in, movement_events);
        end

        bad_chans_fn = get_bad_chan_markings_fn(subject_name, options.tag_bad_channels.export_dir);
        if exist(bad_chans_fn, 'file')
            % TODO: load content of .mat and set channel flag
            % TODO: save data -> see process_nst_tag_bad_channels
        end
    end
    %% Manage bad channel markins
    % TODO: update channel flags
    % 
    files_in{ifile} = file_in;
end

end


function redone=load_anatomy(options,i_subject)

    subject=options.import.subject{i_subject};
    [sSubject, iSubject] = bst_get('Subject',subject.name, 1);

    
    redone=options.import.redo && length(sSubject.Surface) >= 1;
    if redone  || ~size(sSubject.Surface,2)
        fiducial_file=fullfile(subject.mri_folder,'fiducials.mat');
        if exist(fiducial_file,'file')
            fiducials=load(fiducial_file);
            bst_process('CallProcess', 'process_import_anatomy', [], [], ... 
                                   'subjectname', subject.name, ...
                                   'mrifile',     {subject.mri_folder, options.import.mri_folder_type}, ...
                                   'nvertices',   options.import.nvertices, ...
                                   'nas', fiducials.NAS, ...
                                   'lpa', fiducials.LPA, ...
                                   'rpa', fiducials.RPA, ...
                                   'ac' , fiducials.AC,  ...
                                   'pc' , fiducials.PC,  ...
                                   'ih' , fiducials.IH,  ...
                                   'aseg',options.import.aseg );
        else
            bst_process('CallProcess', 'process_import_anatomy', [], [], ... 
                                   'subjectname', subject.name, ...
                                   'mrifile',     {subject.mri_folder, options.import.mri_folder_type}, ...
                                   'nvertices',   options.import.nvertices, ...
                                   'aseg',options.import.aseg );

        end
        
        % Reload subject
        [sSubject, iSubject] = bst_get('Subject',subject.name, 1);

        % Set the default cortex surface
        iCortex = find(strcmp({sSubject.Surface.Comment}, options.head_model.surface));
        db_surface_default(iSubject, 'Cortex', iCortex);
        panel_protocols('RepaintTree');
        
        % Remesh the head surface to have an homogenous mapping of the
        % surface 
        iHead = find( contains( {sSubject.Surface.Comment},'head') );

        HeadSurfaceFile=sSubject.Surface(iHead).FileName;
        [SurfaceMat, SurfaceFile] = in_tess_bst(HeadSurfaceFile);

        tess_remesh(sSubject.Surface(iHead).FileName, size(SurfaceMat.Vertices,1)) 
        
        
        
        % Compute Voroi Interpolator 
        % To fix : It doesn't work if the segmentation file is already
        % imported
        bst_process('CallProcess', 'process_nst_compute_voronoi', [], [], ...
                    'subjectname', subject.name);
        
        
        % Process: Import Segentated MRI
        segmentation_file=fullfile(subject.mri_folder,'Segmentation','segmentation_5tissues.nii');
        if exist(segmentation_file,'file')
            bst_process('CallProcess', 'process_import_mri', [], [], ...
                'subjectname', subject.name, ...
                'mrifile',     {segmentation_file, 'ALL'});
            % Rename and fix segmentation file value 
            
           [sSubject, iSubject] = bst_get('Subject',subject.name, 1);
           iSegmentation = find( contains( {sSubject.Anatomy.Comment},'segmentation') );
           fix_segmentation(sSubject, iSubject,iSegmentation);
        end
        
        
    end
end

function [file_in, redone] = import_nirs_file(options,i_subject)
subject=options.import.subject{i_subject};

%% Import data
nirs_fn = subject.nirs_fn;

if isempty(subject.name)
   [root, subject_name, ext] = fileparts(nirs_fn);
else
    subject_name = subject.name;
end

condition = ['origin' get_ppl_tag()];
[file_in, redone] = nst_run_bst_proc([subject_name '/' condition '/Raw'], options.import.redo, ...
                                       'process_import_data_time', [], [], ...
                                       'subjectname',  subject_name, ...
                                       'condition',    condition, ...
                                       'datafile',     {nirs_fn, 'NIRS-BRS'}, ...
                                       'timewindow',   [], ...
                                       'split',        0, ...
                                       'ignoreshort',  1, ...
                                       'channelalign', 1, ...
                                       'usectfcomp',   0, ...
                                       'usessp',       0, ...
                                       'freq',         [], ...
                                       'baseline',     []);
%% If available import additional headpoints, and refine the registration
if redone && ~isempty(subject.additional_headpoints)  && exist(subject.additional_headpoints)
    % Import head points
    bst_process('CallProcess', 'process_headpoints_add', file_in, [], ...
                'channelfile', {subject.additional_headpoints, 'ASCII_NXYZ'}, ...
                'fixunits',    0.1, ...
                'vox2ras',     1);

    %  Refine registration
    bst_process('CallProcess', 'process_headpoints_refine', file_in, []);
    
end 

%% Manage movement event markings TODO
if redone
    evt_formats = bst_get('FileFilters', 'events');
    evt_format = evt_formats(strcmp('BST', evt_formats(:,3)), :);

    moco_fn = get_moco_markings_fn(subject_name, options.moco.export_dir);
    if exist(moco_fn, 'file')
        % Load event from pre-saved file
        % TODO: test
        sFile_in = load(file_fullpath(file_in));
        [sFile_in, events] = import_events(sFile_in, [], moco_fn, evt_format);   
    else
        % Create empty event group
        movement_events = db_template('event');
        movement_events.label = 'movement_artefacts';
        sFile_in = bst_process('GetInputStruct', file_in);
        process_nst_import_csv_events('import_events', [], sFile_in, movement_events);
    end

    bad_chans_fn = get_bad_chan_markings_fn(subject_name, options.tag_bad_channels.export_dir);
    if exist(bad_chans_fn, 'file')
        % TODO: load content of .mat and set channel flag
        % TODO: save data -> see process_nst_tag_bad_channels
    end
    
end
%% Manage bad channel markins
% TODO: update channel flags
% 

end

function [files_in, redone_imports] = import_subjects(options)
nb_subjects = length(options.import.subject);

files_in = cell(1,nb_subjects);
redone_imports = zeros(1,nb_subjects);

for i=1:nb_subjects
    sSubject = bst_get('Subject', options.import.subject{i}.name, 1);
    if isempty(sSubject)
        [sSubject, iSubject] = db_add_subject(options.import.subject{i}.name, [], 0, 0);
    end

    redone_imports(i)=load_anatomy(options,i);
    options.import.redo = options.import.redo || redone_imports(i);

    db_save();

    [imported_files, redone] = import_nirs_file(options,i);
    files_in{i} = imported_files;
    redone_imports(i)=redone;
end

end



function markings_fn = get_moco_markings_fn(subject_name, export_dir)
markings_fn = '';
if ~isempty(export_dir)
    assert(exist(export_dir, 'dir')~=0);
    markings_fn = fullfile(export_dir, [subject_name '_motion_events.mat']);
end
end

function markings_fn = get_bad_chan_markings_fn(subject_name, export_dir)
markings_fn = '';
if ~isempty(export_dir)
    assert(exist(export_dir, 'dir')~=0);
    markings_fn = fullfile(export_dir, [subject_name '_bad_channels.mat']);
end
end

%% Helper functions

function folder = create_dir(folder)
% Create folder if does not exist. 
% Check that folder is not a subfolder of nirstorm sources (encourage good practice
% not to store data in source code folders)

if ~isempty(folder)
    if exist(fullfile(folder, 'nst_install.m'), 'file') || ...
            exist(fullfile(folder, '..', 'nst_install.m'), 'file') || ...
            exist(fullfile(folder, '..', '..', 'nst_install.m'), 'file')
        warning('Processing folder should not be part of nirstorm source folders (%s)', folder);
    end

    if ~exist(folder, 'dir')
        mkdir(folder);
    end
end

end

function fix_segmentation(sSubject, iSubject,iSegmentation)
% Fix the value of the segmentation MRI file.
% The file should contains value between 0 and 5 instead of 0 and 255 
% This step is needed in order to compute the fluences

segmentation_fn = file_fullpath(sSubject.Anatomy(iSegmentation).FileName);
sSegmentation = load(segmentation_fn);
sSegmentation.Cube(sSegmentation.Cube==51)=1;
sSegmentation.Cube(sSegmentation.Cube==102)=2;
sSegmentation.Cube(sSegmentation.Cube==153)=3;
sSegmentation.Cube(sSegmentation.Cube==204)=4;
sSegmentation.Cube(sSegmentation.Cube==255)=5;

bst_save(segmentation_fn, sSegmentation, 'v7');

db_save();
end



function tag = get_bst_file_tag(fn)

[rr, bfn, ext] = fileparts(fn);
bst_prefixes = {'results_','data_','linkresults_','linkdata_','pdata_','presults_'};
for ipref=1:length(bst_prefixes)
    bfn = replace(bfn, bst_prefixes{ipref}, '');
end

toks = regexp(bfn, '(.*)(?:_\d{6}_\d{4})', 'tokens');
if isempty(toks) || isempty(toks{1})
    tag = bfn;
else
    tag = toks{1}{1};
end

end

function flag = function_exists(func_name)

flag = 1;
try
    eval(func_name);
catch ME
    if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
        flag = 0;
    end
end

end

function sfn = protect_fn_str(s)
sfn = strrep(s, ' | ', '--');
sfn = strrep(s, ' : ', '--');
sfn = strrep(s, ' :', '--');
sfn = strrep(s, ' ', '_');
end


function plot_stat(sFile_ttest, fig_fn, options, show_colbar, save_colbar, sSubjectDefault)

% TODO: set colormap
% TODO: set t-stat cmap boundaries

hFigSurfData = view_surface_data(sSubjectDefault.Surface(sSubjectDefault.iCortex).FileName, ...
                                 sFile_ttest, 'NIRS', 'NewFigure');
StatThreshOptions = bst_get('StatThreshOptions');
StatThreshOptions.pThreshold = options.GLM_group.contrast_tstat.plot.pvalue_threshold;
StatThreshOptions.Correction = options.GLM_group.contrast_tstat.plot.pvalue_mcc_method;
bst_set('StatThreshOptions', StatThreshOptions);

ColormapInfo = getappdata(hFigSurfData, 'Colormap');
ColormapType = ColormapInfo.Type;
switch ColormapType
    case 'stat2'
        bst_colormaps('SetColormapName', 'stat2', 'cmap_ovun');
    case 'stat1'
        bst_colormaps('SetColormapName', 'stat1', 'cmap_ovun');
end
if ~isempty(options.fig_background)
    bst_figures('SetBackgroundColor', hFigSurfData, options.fig_background);
end
if ~show_colbar
    bst_colormaps('SetDisplayColorbar', 'stat2', 0);
end
view(options.fig_cortex_view);
zoom(hFigSurfData, options.fig_cortex_zoom);

iSurf = getappdata(hFigSurfData, 'iSurface');

panel_surface('SetSurfaceSmooth',       hFigSurfData, iSurf, options.fig_cortex_smooth, 0);
% panel_surface('SetSurfaceTransparency', hFigSurfData, iSurf, SurfAlpha);
panel_surface('SetShowSulci',           hFigSurfData, iSurf, options.fig_cortex_show_sulci);
%panel_surface('SetDataThreshold',       hFig, iSurf, DataThreshold);
%panel_surface('SetSizeThreshold',       hFig, iSurf, SizeThreshold);

nst_save_figure(fig_fn, options, hFigSurfData);

if save_colbar
    % Save colorbar
    [root, bfn, ext] = fileparts(fig_fn);
    colbar_fig_fn = fullfile(root, [bfn '_colobar' ext]);
    bst_colormaps('SetDisplayColorbar', 'stat2', 1);
    hColorbar = findobj(GlobalData.CurrentFigure.Last, '-depth', 1, 'Tag', 'Colorbar');
    set(hColorbar, 'XColor', [0 0 0]);
    set(hColorbar, 'YColor', [0 0 0]);
    options_colbar = options;
    options_colbar.export_fig_dpi = 500;
    nst_save_figure(colbar_fig_fn, options_colbar, hColorbar);
end
close(hFigSurfData);

end
