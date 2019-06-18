function varargout = nst_ppl_1st_level_channel_V1(action, options, arg1, arg2)
%NST_PPL_1ST_LEVEL_CHANNEL_V1    
% Run a full channel-space intra-subject pipeline
% starting from raw NIRS data up to GLM contrasts (effect and t-stat maps).
% 
% DEFAULT_OPTIONS = NST_PPL_1ST_LEVEL_CHANNEL_V1('setup', PROTOCOL_NAME)
%     Return default options and set current protocol to PROTOCOL_NAME
%
% [GLM_CON_FILES, REDONE_GLM, PREPROC_FOLDER, GLM_FOLDER] = 
%         NST_PPL_1ST_LEVEL_CHANNEL_V1(FILE_RAW, GLM_EVENTS_ORDER, 
%                                      GLM_CONTRASTS, OPTIONS, FORCE_REDO=0)
%   
%     Apply pipeline to given data FILE_RAW.
%     ASSUME: all subjects in a given protocol have the same template
%             anatomy. TODO: check that and report error if not OK.
%
%     List of steps:
%        1) Export manual inputs: TODO
%           - movement artefacts
%           - bad channels
%        1) Head model for all pairs
%            If needed, clone FILE_RAW and precompute a head model for
%            all possible pairs, which can be used for other subjects.
%            Create a dummy subject called "full_head_model...".
%        1) Motion correction
%            ASSUME: event group "movement_artefacts" exists in FILE_RAW and
%                    has been filled before calling this function.
%                    See script nst_import_data.m TODO
%        2) Resampling:
%            TODO: check interpolation errors when there are spikes 
%        3) Detect bad channels
%        4) Convert to delta optical density
%        5) High pass filter -- TODO: bandpass filter
%        6) Compute head model (from "full_head_model...")
%        7) Project on the cortical surface
%        8) Run GLM:
%            - build design matrix from stimulation events GLM_EVENTS_ORDER
%            - fit GLM
%            - compute contrasts defined in GLM_CONTRASTS
%
% TODO:
% - options documentation
% - check anatomical consistency of head model -> within subhead model
%   process already?
% - export manual inputs
% - importation of manual inputs:
%   NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1('import_manual_markings', PROTOCOL_NAME)
% - wiki page
% - utest
% 
options.PIPELINE_TAG = get_ppl_tag();


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

protocol_info = bst_get('ProtocolInfo');
if protocol_info.UseDefaultAnat~=1
    error('Protocol should use default anatomy for all subjects');
end

create_dir(options.fig_dir);
create_dir(options.moco.export_dir);
create_dir(options.tag_bad_channels.export_dir);
create_dir(options.GLM_group.rois_summary.csv_export_output_dir);

file_raw1 = nst_get_bst_func_files(groups(1).subject_names{1}, ['origin' get_ppl_tag()], 'Raw');
if isempty(file_raw1)
    error(sprintf('Cannot find "origin/Raw" data for subject "%s". Consider using nst_ppl_surface_template_V1(''import'',...).',  ...
                  groups(1).subject_names{1})); %#ok<SPERR>
end

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
        file_raw = nst_get_bst_func_files(subject_name, ['origin' get_ppl_tag()], 'Raw');
        if isempty(file_raw)
            error(sprintf('Cannot find "origin/Raw" data for subject "%s". Consider using nst_ppl_surface_template_V1(''import'',...).',  ...
                          subject_name)); %#ok<SPERR>
        end
        
        % Run preprocessings
        [sFiles_preprocessed{1}, redone_preprocs, preproc_folder_raw, preproc_folder_hb] = preprocs(file_raw,options, force_redo);        
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

    %% Group-level analysis
    if length(subject_names) <= 2
       if ~isempty(group_label)
           group_msg = [' of ' group_label];
       else
           group_msg = '';
       end
       warning('Not enough data for group analysis%s.\n', group_msg);
       return;
    end
    sSubjectDefault = bst_get('Subject', 0); %TODO: use this also for 1st level
    
    hb_types = process_nst_cortical_projection('get_hb_types');
    nb_hb_types=length(hb_types);
    stacking_types = process_nst_concat_matrices('get_stacking_types');
    contrasts = options.GLM_1st_level.contrasts;
    if options.GLM_group.do
       
        if ~isempty(group_label) 
            group_condition_name = [group_label '_'];
        else
            group_condition_name = '';
        end
        group_condition_name = [group_condition_name 'GLM' get_ppl_tag()];
        group_condition_names{igroup} = group_condition_name;
        for ihb=1:nb_hb_types
            for icon=1:nb_contrasts
                item_comment = ['Group analysis/' group_condition_name '/' hb_types{ihb} ' | con_t-+ ' contrasts(icon).label];
                [sFile_GLM_gp_ttest, redone] = nst_run_bst_proc(item_comment, redo_group, ...
                                                                'process_nst_glm_group_ttest', all_sFiles_con{igroup}(ihb, icon, :), [], ...
                                                                'tail', 'two');
                fig_bfn = sprintf('group_%s_%s_tmap_mcc_%s_pv_thresh_%s_%s.png', ...
                              group_label, hb_types{ihb}, options.GLM_group.contrast_tstat.plot.pvalue_mcc_method,...
                              nst_format_pval(options.GLM_group.contrast_tstat.plot.pvalue_threshold), ...
                              contrasts(icon).label);
                fig_fn = fullfile(options.fig_dir, fig_bfn);
                if ~isempty(options.fig_dir) && options.make_figs && ...
                    options.GLM_group.contrast_tstat.plot.do && ...
                    (redone || options.GLM_group.contrast_tstat.plot.redo || ~exist(fig_fn, 'file'))   
                    plot_stat(sFile_GLM_gp_ttest, fig_fn, options, 0, 1, sSubjectDefault);
                end

                if options.GLM_group.rois_summary.do
                    
                    if igroup==1 && ihb==1 && icon==1
                        sFile_gp_masks = cell(nb_groups, nb_hb_types, nb_contrasts);
                    end
                    
                    [sFile_gp_mask, redone] = nst_run_bst_proc(['Group analysis/' group_condition_name '/' hb_types{ihb} ' | con_t-+ ' contrasts(icon).label ' | mask'], ...
                                                                redone | options.GLM_group.rois_summary.redo, ...
                                                                'process_nst_glm_contrast_mask', sFile_GLM_gp_ttest, [], ...
                                                                'do_atlas_inter', 1, ...
                                                                'min_atlas_roi_size', 3, ...
                                                                'atlas', options.GLM_group.rois_summary.atlas);
                     sFile_gp_masks{igroup, ihb, icon} = sFile_gp_mask;
                end   
            end % loop over contrasts
        end % loop over hb types
    end % If do group level    
end % Loop over groups

end


function [sFile_dHb_filtered, redone, preproc_folder_raw, preproc_folder_hb] = preprocs(sFile_raw, options, force_redo)

if nargin < 4
    force_redo = 0;
end

preproc_folder_raw = sprintf('preprocessing_raw_%s/', options.PIPELINE_TAG);

% Compute Scalp coupling index
nst_run_bst_proc([preproc_folder_raw 'SCI'], force_redo | options.sci.redo, 'process_nst_sci', sFile_raw);

% TODO: export motion correction tagging to external file
% sRaw = load(file_fullpath(sFile_raw));
% sExport.Events = sRaw.Events(strcmp({sRaw.Events.label}, 'movement_artefacts'));
% export_events(sExport, [], moco_export_fn);

% TODO: export bad channel tagging information

% Motion correction
redo_parent = force_redo | options.moco.redo;
[sFileMoco, redo_parent] = nst_run_bst_proc([preproc_folder_raw 'Motion-corrected'], redo_parent, ...
                             'process_nst_motion_correction', sFile_raw, [], ...
                             'option_event_name', 'movement_artefacts');

                         
% Resample to 5Hz (save some space)
redo_parent = redo_parent | options.resample.redo;
[sFileMocoResampled, redo_parent] = nst_run_bst_proc([preproc_folder_raw 'Motion-corrected | Resampled'], redo_parent, ...
                                                      'process_resample', sFileMoco, [], ...
                                                      'freq', options.resample.freq, ...
                                                      'read_all', 1);
% Process: Detect bad channels
% This one is done in-place -> not tracked to handle all do/redo scenarios
if redo_parent
    bst_process('CallProcess', 'process_nst_detect_bad', sFileMocoResampled, [], ...
                'option_remove_negative', 1, ...
                'option_invalidate_paired_channels', 1, ...
                'option_max_sat_prop', options.tag_bad_channels.max_prop_sat_ceil, ...
                'option_min_sat_prop', options.tag_bad_channels.max_prop_sat_floor);
end

% Convert to Hb
preproc_folder_hb = sprintf('preprocessing_Hb_%s/', options.PIPELINE_TAG);
redo_parent = redo_parent | options.MBLL.redo;
[sFile_dHb, redo_parent] = nst_run_bst_proc([preproc_folder_hb 'dHb'], redo_parent, 'process_nst_mbll', sFileMocoResampled, [], ...
                                            'option_age',             options.MBLL.age, ...
                                            'option_pvf',             options.MBLL.pvf, ...
                                            'option_do_plp_corr',     1, ...
                                            'option_dpf_method',      2, ...  % DUNCAN1996
                                            'option_baseline_method', options.MBLL.baseline_method, ...  % median
                                            'timewindow',             []);
              
% Band pass filter
redo_parent = redo_parent | options.high_pass_filter.redo;
[sFile_dHb_filtered, redo_parent] = nst_run_bst_proc([preproc_folder_hb 'Hb | filtered'],  redo_parent, 'process_bandpass', sFile_dHb, [], ...
                                                     'highpass', {options.high_pass_filter.low_cutoff, ''}, ...
                                                     'lowpass', {0, ''}, ...
                                                     'attenuation', 'relax', ...
                                                     'mirror', 0, ...
                                                     'sensortypes', 'NIRS');
                                  
redone = redo_parent;
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
                                                       'fitting',        1, ...  % OLS 
                                                       'statistical_processing', 1, ... % precoloring
                                                       'hpf_low_cutoff', options.high_pass_filter.low_cutoff, ...
                                                       'stim_events',    strjoin(stim_events, ', '), ...
                                                       'hrf_model',      1, ...  % CANONICAL
                                                       'trend',          1, ...
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
options.import.subject{1}.name='';
options.import.subject{1}.nirs_fn='';
options.import.subject{1}.additional_headpoints='';
options.sci.redo = 0;


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

options.MBLL.redo = 0;
options.MBLL.age = 60;
options.MBLL.pvf = 50;
options.MBLL.baseline_method = 2; % 1:mean, 2:median


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


function options = setup(protocol_name)

if nargin >= 1
    sProtocol = bst_get('Protocol', protocol_name);
    if isempty(sProtocol)
        error('Protocol %s not found', protocol_name);
    end
    gui_brainstorm('SetCurrentProtocol', sProtocol);
end

options.sci.redo = 0;

options.moco.redo = 0;
options.moco.export_dir = fullfile('.', 'moco_marking');

options.high_pass_filter.redo = 0;
options.high_pass_filter.low_cutoff = 0.01; %Hz

options.tag_bad_chans.redo = 0;
options.tag_bad_channels.max_prop_sat_ceil = 1; % no tagging
options.tag_bad_channels.max_prop_sat_floor = 1; % no tagging

options.GLM.redo = 0;
options.GLM.trim_start = 0; % sec
options.GLM.save_fit = 0;

options.contrast.redo = 0;
options.contrast_plot.redo = 0;

options.make_figs = 1;
options.fig_print_method = 'saveas'; % 'saveas', 'export_fig'
options.fig_dir = fullfile('.', 'figs');
options.export_fig_dpi = 90;

options.clean_preprocessings = 0;

% Oblique view from the top
% options.plot_3d_view_az = ;
% options.plot_3d_view_el = ;

if ~function_exists('export_fig')
    disp('"export_fig" not found. Can be installed from "https://github.com/altmany/export_fig"');
    return
end

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

    db_save();

    [imported_files, redone] = import_nirs_file(options,i);
    files_in{i} = imported_files;
    redone_imports(i)=redone;
end

end

function hplots = plot_paradigm(paradigm, y_min, y_max, cond_colors)

hold on;
hplots = zeros(length(paradigm), 1);
for icond=1:length(paradigm)
    for itrial=1:size(paradigm(icond).times, 2)
        trial_start = paradigm(icond).times(1, itrial);
        trial_end = paradigm(icond).times(2, itrial);
        hplots(icond) = area([trial_start trial_end],  [y_max y_max], ...
                             y_min, 'FaceColor', cond_colors(icond, :), ...
                             'EdgeColor', [1, 1, 1]);                  
    end
end

end


%% Helper functions
function ptag = get_ppl_tag()
ptag = '_nspc_V1';
end

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
