function varargout = nst_ppl_surface_template_V1(action, options, arg1, arg2)
%NST_PPL_SURFACE_TEMPLATE_V1
% Manage a full template- and surface-based pipeline starting from raw NIRS data
% up to GLM group analysis (if enough subjects).
% Can keep track of user-defined markings outside of brainstorm db 
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
    case 'import'
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

% Set default cortical surface
sSubject = bst_get('Subject', 0);
prev_iCortex = sSubject.iCortex;
%TODO: only if necessary
iCortex = find(strcmp({sSubject.Surface.Comment}, options.head_model.surface));
db_surface_default(0, 'Cortex', iCortex);
panel_protocols('RepaintTree');

create_dir(options.fig_dir);
create_dir(options.moco.export_dir);
create_dir(options.tag_bad_channels.export_dir);
create_dir(options.GLM_group.rois_summary.csv_export_output_dir);

% Get head model precomputed for all optode pairs
% (precompute it by cloning given data if needed)
file_raw1 = nst_get_bst_func_files(groups(1).subject_names{1}, ['origin' get_ppl_tag()], 'Raw');
if isempty(file_raw1)
    error(sprintf('Cannot find "origin/Raw" data for subject "%s". Consider using nst_ppl_surface_template_V1(''import'',...).',  ...
                  groups(1).subject_names{1})); %#ok<SPERR>
end
[sFile_raw_fhm, fm_subject, fhm_redone] = get_sFile_for_full_head_model(file_raw1, options, force_redo);

%% Within-subject analyses
nb_groups = length(groups);
all_sFile_table_zscores = cell(1, nb_groups);
group_condition_names = cell(1, nb_groups);
any_rois_summary_redone = 0;
for igroup=1:nb_groups
    redo_group = options.GLM_group.redo;
    subject_names = groups(igroup).subject_names;
    group_label = groups(igroup).label;
    for isubject=1:length(subject_names)

        subject_name = subject_names{isubject};

        % Check if given subject is dummy one created to hold full head model    
        if strcmp(subject_name, fm_subject)
            warning('Ignoring dummy subject "%s" used to store head model for all pairs', fm_subject);
            continue;
        end

        file_raw = nst_get_bst_func_files(subject_name, ['origin' get_ppl_tag()], 'Raw');
        if isempty(file_raw)
            error(sprintf('Cannot find "origin/Raw" data for subject "%s". Consider using nst_ppl_surface_template_V1(''import'',...).',  ...
                          subject_name)); %#ok<SPERR>
        end
        
        % Run preprocessings
        [sFiles_preprocessed, hb_types, redone_preprocs, preproc_folder] = preprocs(file_raw, sFile_raw_fhm, options, force_redo|fhm_redone);

        % Run 1st level GLM
        [sFiles_GLM, sFiles_con, redone_any_contrast, glm_folder] = glm_1st_level(sFiles_preprocessed, options, ...
                                                                                  redone_preprocs | force_redo);

        if isubject==1
            all_sFiles_con = cell(size(sFiles_con, 1), size(sFiles_con, 2), length(subject_names));
        end
        all_sFiles_con(:,:,isubject) = sFiles_con;

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
    all_sFiles_subj_zmat = cell(size(all_sFiles_con, 1), size(all_sFiles_con, 2));
    all_sFiles_subj_zmat_con_concat = cell(1, size(all_sFiles_con, 1));

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
        for ihb=1:size(all_sFiles_con, 1)
            for icon=1:size(all_sFiles_con, 2)
                % TODO: use pipeline tag to name condition folder
                item_comment = ['Group analysis/' group_condition_name '/' hb_types{ihb} ' | con_t-+ ' contrasts(icon).label];
                [sFile_GLM_gp_ttest, redone] = nst_run_bst_proc(item_comment, redo_group, ...
                                                                'process_nst_glm_group_ttest', all_sFiles_con(ihb, icon, :), [], ...
                                                                'tail', 'two');
                fig_bfn = sprintf('group_%s_%s_tmap_mcc_%s_pv_thresh_%s_%s.png', ...
                              group_label, hb_types{ihb}, options.GLM_group.contrast_tstat.plot.pvalue_mcc_method,...
                              nst_format_pval(options.GLM_group.contrast_tstat.plot.pvalue_threshold), ...
                              contrasts(icon).label);
                fig_fn = fullfile(options.fig_dir, fig_bfn);
                if options.make_figs && options.GLM_group.contrast_tstat.plot.do && ...
                    (redone || options.GLM_group.contrast_tstat.plot.redo || ~exist(fig_fn, 'file'))
                    hFigSurfData = view_surface_data(sSubjectDefault.Surface(sSubjectDefault.iCortex).FileName, ...
                                                     sFile_GLM_gp_ttest, 'NIRS', 'NewFigure');
                    StatThreshOptions = bst_get('StatThreshOptions');
                    StatThreshOptions.pThreshold = options.GLM_group.contrast_tstat.plot.pvalue_threshold;
                    StatThreshOptions.Correction = options.GLM_group.contrast_tstat.plot.pvalue_mcc_method;
                    bst_set('StatThreshOptions', StatThreshOptions);

                    % TODO: set surface smoothing
                    % TODO: set better colormap that does not span values betwn -3 and 3
                    % TODO: factorize map plotting
                    if ~isempty(options.fig_background)
                        bst_figures('SetBackgroundColor', hFigSurfData, options.fig_background);
                    end
                    bst_colormaps('SetDisplayColorbar', 'stat2', 0);
                    view(options.fig_cortex_view);
                    zoom(hFigSurfData, options.fig_cortex_zoom);
                    nst_save_figure(fig_fn, options, hFigSurfData);

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

                    close(hFigSurfData);
                end

                if options.GLM_group.rois_summary.do
                    [sFile_gp_mask, redone] = nst_run_bst_proc(['Group analysis/' group_condition_name '/' hb_types{ihb} ' | con_t-+ ' contrasts(icon).label ' | mask'], ...
                                                                redone | options.GLM_group.rois_summary.redo, ...
                                                                'process_nst_glm_contrast_mask', sFile_GLM_gp_ttest, [], ...
                                                                'do_atlas_inter', 1, ...
                                                                'min_atlas_roi_size', 3, ...
                                                                'atlas', options.GLM_group.rois_summary.atlas);
                     % TODO: switch between matrix output (atlas-based) or map output (full mask)                                       
                     [sFile_subj_zmat, redone] = nst_run_bst_proc(['Group analysis/' group_condition_name '/' hb_types{ihb} ' | con ' contrasts(icon).label ' | masked z-scores'], ...
                                                                   redone | options.GLM_group.rois_summary.redo, ...
                                                                   'process_nst_glm_group_subjs_zmat', ...
                                                                   all_sFiles_con(ihb, icon, :), sFile_gp_mask);
                     all_sFiles_subj_zmat{ihb, icon} = sFile_subj_zmat;
                end
                
                if redone
                    % Set contrast name as prefix for each ROI column                   
                    bst_process('CallProcess', 'process_nst_prefix_matrix', ...
                                sFile_subj_zmat, [], 'col_prefixes', [contrasts(icon).label '_']);
                end
                
            end
            if options.GLM_group.rois_summary.do
                % Concatenate across contrasts
                [sFile_concat, redone] = nst_run_bst_proc(['Group analysis/' group_condition_name '/' hb_types{ihb} ' masked z-scores'], ...
                                                         redone | options.GLM_group.rois_summary.redo, ...
                                                         'process_nst_concat_matrices', ...
                                                          all_sFiles_subj_zmat(ihb, :), [], ...
                                                         'stacking_type', stacking_types.column);
                all_sFiles_subj_zmat_con_concat{ihb} = sFile_concat;
                
                if redone
                    % Set Hb type as prefix for all columns. Add user-defined col prefix
                    if ~isempty(options.GLM_group.rois_summary.matrix_col_prefix)
                        col_prefix = [options.GLM_group.rois_summary.matrix_col_prefix '_'];
                    else
                        col_prefix = '';
                    end
                    bst_process('CallProcess', 'process_nst_prefix_matrix', ...
                                sFile_concat, [], 'col_prefixes', [col_prefix hb_types{ihb} '_']);
                end
            end
        end

        if options.GLM_group.rois_summary.do

            [sFile_table_zscores, redone] = nst_run_bst_proc(['Group analysis/' group_condition_name '/all masked zscores'], ...
                                                         redone | options.GLM_group.rois_summary.redo, ...
                                                         'process_nst_concat_matrices', ...
                                                          all_sFiles_subj_zmat_con_concat, [], ...
                                                         'stacking_type', stacking_types.column);

            all_sFile_table_zscores{igroup} = sFile_table_zscores;
            any_rois_summary_redone = any_rois_summary_redone | redone;
            

            
            if redone && isempty(options.GLM_group.rois_summary.stack_groups)
                % Group results will not be stacked so save each group data separately
                if isempty(group_label)
                    group_prefix = '';
                else
                    group_prefix = [group_label '_'];
                end
                csv_fn = fullfile(options.GLM_group.rois_summary.csv_export_output_dir, ...
                                  [group_prefix 'z-scores.csv']);
                bst_process('CallProcess', 'process_nst_save_matrix_csv', ...
                            sFile_table_zscores, [], ...
                            'ignore_rows_all_zeros', 0, 'ignore_cols_all_zeros', 0, ...
                            'csv_file', {csv_fn, 'ASCII-CSV'});
            end
        end
    end    
end

if options.GLM_group.do && options.GLM_group.rois_summary.do && ...
        ~isempty(options.GLM_group.rois_summary.stack_groups)
    
    % Concatenate all group results and save as CSV
    
    % TODO: move checks to global option checks
    assert(size(options.GLM_group.rois_summary.stack_groups, 1) == 1);
    assert(length(unique(options.GLM_group.rois_summary.stack_groups)) == length(options.GLM_group.rois_summary.stack_groups));
    assert(isempty(setdiff(options.GLM_group.rois_summary.stack_groups, 1:nb_groups)));
    
    [varying, common_pref, common_suf] = str_remove_common(group_condition_names);
    stacked_group_cond_name = [common_pref strjoin(varying, '_') common_suf];
    [sFile_table_zscores, redone] = nst_run_bst_proc(['Group analysis/' stacked_group_cond_name '/all groups masked zscores'], ...
                                                        any_rois_summary_redone | options.GLM_group.rois_summary.redo, ...
                                                        'process_nst_concat_matrices', ...
                                                        all_sFile_table_zscores(options.GLM_group.rois_summary.stack_groups), [], ...
                                                        'stacking_type', stacking_types.row);

    [varying_label, common_prefix, common_suffix] = str_remove_common({groups(options.GLM_group.rois_summary.stack_groups).label}, 1);
    varying_label(cellfun(@isempty, varying_label)) = {''};
    
    csv_fn = fullfile(options.GLM_group.rois_summary.csv_export_output_dir, ...
                      [common_prefix strjoin(varying_label, '_') common_suffix '_z-scores.csv']);
    nst_run_bst_proc({}, redone | options.GLM_group.rois_summary.redo, ...
        'process_nst_save_matrix_csv', ...
        sFile_table_zscores, [], ...
        'ignore_rows_all_zeros', 0, 'ignore_cols_all_zeros', 1, ...
        'csv_file', {csv_fn, 'ASCII-CSV'});
end


%% Finalize
if prev_iCortex ~= iCortex
   % Set default cortical surface to original one
   db_surface_default(0, 'Cortex', prev_iCortex);
   panel_protocols('RepaintTree');
end

% if nargout >= 1
%     varargout{1} = sFiles_con;
% end
% 
% if nargout >= 2
%     varargout{2} = redone_any_contrast;
% end
% 
% if nargout >= 3
%     varargout{3} = sFiles_GLM;
% end
% 
% if nargout >= 4 
%     varargout{4} = glm_folder;
% end
% 
% if nargout >= 5 && ~options.clean_preprocessings
%     varargout{5} = preproc_folder;
% end

end


function [sFilesHbProj, hb_types, redone, preproc_folder] = preprocs(sFile_raw, sFile_raw_full_head_model, options, force_redo)

if nargin < 4
    force_redo = 0;
end

preproc_folder = sprintf('preprocessing%s/', get_ppl_tag());

% Compute Scalp coupling index
nst_run_bst_proc([preproc_folder 'SCI'], force_redo | options.sci.redo, 'process_nst_sci', sFile_raw);

% TODO: export motion correction tagging to external file
% sRaw = load(file_fullpath(sFile_raw));
% sExport.Events = sRaw.Events(strcmp({sRaw.Events.label}, 'movement_artefacts'));
% export_events(sExport, [], moco_export_fn);

% TODO: export bad channel tagging information

% Motion correction
redo_parent = force_redo | options.moco.redo;
[sFileMoco, redo_parent] = nst_run_bst_proc([preproc_folder 'Motion-corrected'], redo_parent, ...
                             'process_nst_motion_correction', sFile_raw, [], ...
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
[dummy_out, redo_parent] = nst_run_bst_proc([preproc_folder 'head model'], redo_parent, 'process_nst_sub_headmodel', ...
                                            sFile_dOD_filtered, sFile_raw_full_head_model);                                  
                                  
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
            if options.make_figs && options.GLM_1st_level.contrast_tstat.plot.do && ...
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

function [file_raw_fm, fm_subject, redone] = get_sFile_for_full_head_model(sfile_raw, options, force_redo)

redone = 0;
fm_subject = ['full_head_model' get_ppl_tag()];
subject_name = fileparts(sfile_raw);
% Check if given subject is dummy one created to hold full head model
if strcmp(subject_name, fm_subject)
    file_raw_fm = sfile_raw;
    if ~force_redo
        return;
    end
else 
    file_raw_fm = nst_get_bst_func_files(fm_subject, ['origin' get_ppl_tag()], 'Raw');
end

head_model_comment = 'head model [all pairs]';

if isempty(file_raw_fm)
    [SubjectName, origin_folder] = bst_fileparts(bst_fileparts(sfile_raw), 1);

    % Lazy way of duplicating data along with channel definition
    % -> export as .nirs, then reimport
    tmp_dir = tempname;
    mkdir(tmp_dir);
    [o1, o2, sInputs] = bst_process('CallProcess', 'process_nst_export_nirs', sfile_raw, [], ...
                                    'outputdir', {tmp_dir, {}});
    nirs_bfn = process_nst_export_nirs('get_output_nirs_bfn', sInputs);
    
    tmp_nirs_fn = fullfile(tmp_dir, nirs_bfn);
    
    sFile_in = bst_process('CallProcess', 'process_import_data_time', [], [], ...
                            'subjectname',  fm_subject, ...
                            'condition',    origin_folder, ...
                            'datafile',     {tmp_nirs_fn, 'NIRS-BRS'}, ...
                            'timewindow',   [], ...
                            'split',        0, ...
                            'ignoreshort',  1, ...
                            'channelalign', 1, ...
                            'usectfcomp',   0, ...
                            'usessp',       0, ...
                            'freq',         [], ...
                            'baseline',     []);
    file_raw_fm = bst_process('CallProcess', 'process_set_comment', sFile_in, [], ...
                            'tag',     'Raw', ...
                            'isindex', 0);
    if isstruct(file_raw_fm)
        file_raw_fm = file_raw_fm.FileName;
    end
    rmdir(tmp_dir, 's');
else
    % TODO: Check that existing surface of existing head model is consistent with
    % options
%     sStudy = bst_get('StudyWithCondition', fileparts(file_raw_fm));
%     parent_head_model_fn = sStudy.HeadModel(1).FileName;
%     parent_head_model = in_bst_headmodel(parent_head_model_fn);
end

% Compute head model for all pairs if needed
[dummy_out, redone] = nst_run_bst_proc(head_model_comment, options.head_model.redo || force_redo, ...
                                       'process_nst_import_head_model', file_raw_fm, [], ...
                                       'use_closest_wl', 1, 'use_all_pairs', 1, ...
                                       'force_median_spread', 0, ...
                                       'normalize_fluence', 1, ...
                                       'smoothing_fwhm', 0);

end

function options = get_options()

options.redo_all = 0;

options.import.redo = 0;

options.head_model.surface = 'cortex_lowres';

options.sci.redo = 0;

options.head_model.redo = 0;

options.moco.redo = 0;
options.moco.export_dir = fullfile('.', 'moco_marking');

options.resample.redo = 0;
options.resample.freq = 5; % Hz

options.dOD.redo = 0;
options.dOD.baseline_def = 0; % 0: mean, 1: median

options.high_pass_filter.redo = 0;
options.high_pass_filter.low_cutoff = 0.01; %Hz

options.tag_bad_channels.redo = 0;
options.tag_bad_channels.max_prop_sat_ceil = 1; % no tagging
options.tag_bad_channels.max_prop_sat_floor = 1; % no tagging
options.tag_bad_channels.export_dir = fullfile('.', 'moco_marking');

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
options.GLM_group.rois_summary.do = 0;
options.GLM_group.rois_summary.atlas = 'MarsAtlas';
options.GLM_group.rois_summary.matrix_col_prefix = '';
options.GLM_group.rois_summary.csv_export_output_dir = 'results';
options.GLM_group.rois_summary.stack_groups = [];

options.make_figs = 1;
options.save_fig_method = 'saveas'; % 'saveas', 'export_fig'
options.export_fig_dpi = 90;
options.fig_dir = fullfile('.', 'figs');
options.fig_background = []; % use default background
options.fig_cortex_view = [89 -24]; % Azimuth and Elevation
                                    % to adjust them manually, right-click on fig 
                                    % then Figure > Matlab controls
                                    % use rotate 3D tool, while moving Az and El
                                    % are displayed in the bottom right of the figure.
options.fig_cortex_zoom = 1;

% Oblique view from the top
% options.plot_3d_view_az = ;
% options.plot_3d_view_el = ;
end

function ptag = get_ppl_tag()
ptag = '__nspst_V1';
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

if exist(fullfile(folder, 'nst_install.m'), 'file') || ...
        exist(fullfile(folder, '..', 'nst_install.m'), 'file') || ...
        exist(fullfile(folder, '..', '..', 'nst_install.m'), 'file')
    warning('Data folder should not be part of nirstorm source folders (%s)', folder);
end

if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
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