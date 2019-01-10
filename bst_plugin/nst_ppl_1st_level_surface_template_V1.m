function varargout = nst_ppl_1st_level_surface_template_V1(arg1, arg2, glm_contrasts, ...
                                                           options, force_redo)
%NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1    
% Run a full template- and surface-based intra-subject pipeline
% starting from raw NIRS data up to GLM contrasts (effect and optional t-stat maps).
% 
% DEFAULT_OPTIONS = NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1('setup')
%     Return default options
%
% DEFAULT_OPTIONS = NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1('setup', PROTOCOL_NAME)
%     Return default options and set current protocol to PROTOCOL_NAME
%
% [GLM_CON_FILES, REDONE_GLM, PREPROC_FOLDER, GLM_FOLDER] = 
%         NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1(FILE_RAW, GLM_EVENTS_ORDER, 
%                                               GLM_CONTRASTS, OPTIONS, FORCE_REDO=0)
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
% - export manual inputs
% - importation of manual inputs:
%   NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1('import_manual_markings', PROTOCOL_NAME)
% - wiki page
% - utest
% 
options.PIPELINE_TAG = '_nspst_V1';

if ischar(arg1) && strcmp(arg1, 'setup')
    if nargin > 1
        assert(ischar(arg2));
        varargout{1} = setup(arg2);
    else
        varargout{1} = setup();
    end
    return;
else
    file_raw = arg1;
    glm_events_order = arg2;
end

if isempty(file_raw)
    error('Empty input file');
end

if nargin < 5
    force_redo = 0;
end

create_dir(options.fig_dir);
create_dir(options.moco.export_dir);

if strcmp(options.save_fig_method, 'export_fig') && ~function_exists('export_fig')
    error('"export_fig" not found. Can be installed from "https://github.com/altmany/export_fig"');
end

% Check if given subject is dummy one created to hold full head model
subject_name = fileparts(fileparts(file_raw));
fm_subject = ['full_head_model_' options.PIPELINE_TAG];
if strcmp(subject_name, fm_subject) % TODO factorize filename logic
    warning('Ignoring dummy subject "%s" used to store head model for all pairs', fm_subject);
    varargout = {}; 
    return;
end

protocol_info = bst_get('ProtocolInfo');
assert(protocol_info.UseDefaultAnat==1);

% Set default cortical surface
sSubject = bst_get('Subject', 0);
prev_iCortex = sSubject.iCortex;
iCortex = find(strcmp({sSubject.Surface.Comment}, options.head_model.surface));
db_surface_default(0, 'Cortex', iCortex);
panel_protocols('RepaintTree');


% Get head model precomputed for all optode pairs
% (precompute it by cloning given data if needed)
[sFile_raw_head_model, fhm_redone] = get_sFile_for_full_head_model(file_raw, options, force_redo);
[sFiles_preprocessed, redone_preprocs, preproc_folder] = preprocs(file_raw, sFile_raw_head_model, options, force_redo|fhm_redone);

[sFiles_GLM, sFiles_con, redone_any_contrast, glm_folder] = glm_1st_level(sFiles_preprocessed, glm_events_order, glm_contrasts, options, redone_preprocs | force_redo);

if options.clean_preprocessings
    full_preproc_folder = fileparts(sFiles_preprocessed{1});
    [sStudy, iStudy] = bst_get('StudyWithCondition', full_preproc_folder);
    db_delete_studies(iStudy);
end

if nargout >= 1
    varargout{1} = sFiles_con;
end

if nargout >= 2
    varargout{2} = redone_any_contrast;
end

if nargout >= 3
    varargout{3} = sFiles_GLM;
end

if nargout >= 4 
    varargout{4} = glm_folder;
end

if nargout >= 5 && ~options.clean_preprocessings
    varargout{5} = preproc_folder;
end

end


function [sFilesHbProj, redone, preproc_folder] = preprocs(sFile_raw, sFile_raw_full_head_model, options, force_redo)

if nargin < 4
    force_redo = 0;
end

preproc_folder = sprintf('preprocessing_%s/', options.PIPELINE_TAG);

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
% This one is done in-place -> not tracked to handle all do/redo scenarios
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
proj_methods = process_nst_cortical_projection('methods');
%TODO: expose as option?
proj_method = proj_methods.Sensitivity_based_interpolation;
% proj_method = proj_methods.MNE;
[sFilesHbProj, redo_parent] = nst_run_bst_proc({[preproc_folder 'dHbO_cortex'], [preproc_folder 'dHbR_cortex']},  redo_parent, ... 
                                                'process_nst_cortical_projection', sFile_dOD_filtered, [], ...
                                                'method', proj_method, ...
                                                'sparse_storage', options.projection.sparse_storage);
redone = redo_parent;
end

function [sFiles_GLM, sFiles_con, redone_any_contrast, glm_folder] = glm_1st_level(sFiles, stim_events, contrasts, options, force_redo)

if nargin < 5
    force_redo = 1;
end
glm_folder = sprintf('GLM_%s/', options.PIPELINE_TAG);

[SubjectName, preprocs_folder] = bst_fileparts(bst_fileparts(sFiles{1}), 1);
sSubject = bst_get('Subject', SubjectName);

redone_any_contrast = 0; % Track if any contrast for any file had to be recomputed
sFiles_GLM = cell(1, length(sFiles));

redo_parent = force_redo | options.GLM.redo;
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
                                                       'trim_start', options.GLM.trim_start, ...
                                                       'save_residuals', 0, ...
                                                       'save_betas',     0, ...
                                                       'save_fit',       0);
end

sFiles_con = cell(length(sFiles_GLM), length(contrasts));
for ifile=1:length(sFiles_GLM)

    for icon=1:length(contrasts)

        % Process: GLM - intra subject contrast
        redo = redone_fit | options.GLM.contrast.redo;
        
        [sFile_GLM_con, redone_con] = nst_run_bst_proc([comment_glm_prefix{ifile} ' | con ' contrasts(icon).label], redo, ...
                                                       'process_nst_glm_contrast', sFiles_GLM{ifile}, [], ...
                                                       'Contrast', contrasts(icon).vector);
        sFiles_con{ifile,icon} = sFile_GLM_con;
        
        % GLM - tmaps
        if options.GLM.contrast_tstat.do
            redo = redone_con | options.GLM.contrast_tstat.redo;
            sFile_GLM_ttest = nst_run_bst_proc([comment_glm_prefix{ifile} ' | con_t-+ ' contrasts(icon).label], redo, ...
                                               'process_nst_glm_contrast_ttest', sFile_GLM_con, [], ...
                                               'tail', 'two');
        
            % Plots
            data_tag = get_bst_file_tag(sFiles{ifile});
            fig_bfn = sprintf('%s_%s_tmap_mcc_%s_pv_thresh_%s_%s.png', ...
                SubjectName, data_tag, options.GLM.contrast_tstat.plot.pvalue_mcc_method,...
                nst_format_pval(options.GLM.contrast_tstat.plot.pvalue_threshold), ...
                contrasts(icon).label);
            fig_fn = protect_fn_str(fullfile(options.fig_dir, fig_bfn ));
            if options.make_figs && options.GLM.contrast_tstat.plot.do && ...
                    (redo || options.GLM.contrast_tstat.plot.redo || ~exist(fig_fn, 'file'))
                hFigSurfData = view_surface_data(sSubject.Surface(sSubject.iCortex).FileName, ...
                    sFile_GLM_ttest, 'NIRS', 'NewFigure');
                StatThreshOptions = bst_get('StatThreshOptions');
                StatThreshOptions.pThreshold = options.GLM.contrast_tstat.plot.pvalue_threshold;
                StatThreshOptions.Correction = options.GLM.contrast_tstat.plot.pvalue_mcc_method;
                %StatThreshOptions.Control    = [1 2 3]; % ???
                bst_set('StatThreshOptions', StatThreshOptions);
            
                % TODO: set surface smoothing
                % TODO: set better colormap that does not span values betwn -3 and 3
                if ~isempty(options.fig_background)
                    bst_figures('SetBackgroundColor', hFigSurfData, options.fig_background);
                end
                bst_colormaps('SetDisplayColorbar', 'stat2', 0);
                view(options.fig_cortex_view); %TODO: expose as option
                zoom(hFigSurfData, options.fig_cortex_zoom);
                nst_save_figure(fig_fn, options, hFigSurfData);
                close(hFigSurfData);
            end
        end
        
        redone_any_contrast = redone_any_contrast | redone_con;
    end
end
end

function [file_raw_fm, redone] = get_sFile_for_full_head_model(sfile_raw, options, force_redo)

fm_subject = ['full_head_model_' options.PIPELINE_TAG];
file_raw_fm = nst_get_bst_func_files(fm_subject, 'origin', 'Raw');

head_model_comment = 'head model [all pairs]';

if isempty(file_raw_fm)
    [SubjectName, origin_folder] = bst_fileparts(bst_fileparts(sfile_raw), 1);

    % Lazy way of duplicating data along with channel definition
    % -> export as .nirs, then reimport
    tmp_dir = tempname;
    mkdir(tmp_dir);
    bst_process('CallProcess', 'process_nst_export_nirs', sfile_raw, [], ...
                'outputdir', {tmp_dir, {}});    
    tmp_nirs_fn = fullfile(tmp_dir, [origin_folder '_Raw.nirs']); %TODO: factorize logic for forging file name
    
    sFile_in = bst_process('CallProcess', 'process_import_data_time', [], [], ...
                            'subjectname',  fm_subject, ...
                            'condition',    'origin', ...
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
                                       'force_median_spread', 1, ...
                                       'normalize_fluence', 0, ...
                                       'smoothing_fwhm', 0);

end

function options = setup(protocol_name)

if nargin >=1
    sProtocol = bst_get('Protocol', protocol_name);
    if isempty(sProtocol)
        error('Protocol %s not found', protocol_name);
    end
    gui_brainstorm('SetCurrentProtocol', sProtocol);
end

options.head_model.surface = 'white_lowres';

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

options.tag_bad_chans.redo = 0;
options.tag_bad_channels.max_prop_sat_ceil = 1; % no tagging
options.tag_bad_channels.max_prop_sat_floor = 1; % no tagging

options.projection.redo = 0;
options.projection.sparse_storage = 0;

options.GLM.redo = 0;
options.GLM.trim_start = 0; % sec

options.GLM.contrast.redo = 0;

options.GLM.contrast_tstat.do = 0; % not active by default -> only beta values are mandatory for group-level analysis
options.GLM.contrast_tstat.redo = 0;
options.GLM.contrast_tstat.plot.do = 0; % not active by default -> can produce a lot of figures
options.GLM.contrast_tstat.plot.redo = 0; 
options.GLM.contrast_tstat.plot.pvalue_threshold = 0.05;
options.GLM.contrast_tstat.plot.pvalue_mcc_method = 'none';

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

options.clean_preprocessings = 0;

% Oblique view from the top
% options.plot_3d_view_az = ;
% options.plot_3d_view_el = ;
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

if ~exist(folder, 'dir')
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