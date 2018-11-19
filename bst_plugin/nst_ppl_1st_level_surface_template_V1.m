function varargout = nst_ppl_1st_level_surface_template_V1(arg1, arg2, glm_contrasts, ...
                                                           options, force_redo)
%NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1    
% Run a full template- and surface-based intra-subject pipeline
% starting from raw NIRS data up to GLM contrasts (effect and t-stat maps).
% 
% DEFAULT_OPTIONS = NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1('setup', PROTOCOL_NAME)
%     Return default options and set current protocol to PROTOCOL_NAME
%
% [GLM_CON_FILES, REDONE] = NST_PPL_1ST_LEVEL_SURFACE_TEMPLATE_V1(FILE_RAW, GLM_EVENTS_ORDER, GLM_CONTRASTS, OPTIONS, FORCE_REDO=0)
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
% - check anatomical consistency of head model -> within subhead model
%   process
% - export manual inputs
% - wiki page
% - utest
% - allow sparse storage for projected signals
%
options.PIPELINE_TAG = '_nspst_V1';

if ischar(arg1) && strcmp(arg1, 'setup')
    assert(ischar(arg2));
    varargout{1} = setup(arg2);
    return;
else
    file_raw = arg1;
    glm_events_order = arg2;
end

if nargin < 5
    force_redo = 0;
end

if ~isempty(file_raw)
    % Get head model precomputed for all optode pairs
    % (precompute it by cloning given data if needed)
    sFile_raw_head_model = get_full_head_model(file_raw, options);
    [sFiles_preprocessed, redone_preprocs] = preprocs(file_raw, sFile_raw_head_model, options);
else
    error('Given raw data is empty');
end

[varargout{1:nargout}] = glm_1st_level(sFiles_preprocessed, glm_events_order, glm_contrasts, options, redone_preprocs | force_redo);

end


function [sFilesHbProj, redo_parent] = preprocs(sFile_raw, sFile_raw_full_head_model, options, force_redo)

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
sFileMoco = nst_run_bst_proc([preproc_folder 'Motion-corrected'], redo_parent, ...
                             'process_nst_motion_correction', sFile_raw, [], ...
                             'option_event_name', 'movement_artefacts');

                         
% Resample to 5Hz (save some space)
redo_parent = redo_parent | options.resample.redo;
sFileMocoResampled = nst_run_bst_proc([preproc_folder 'Motion-corrected | Resampled'], redo_parent, ...
                                      'process_resample', sFileMoco, [], ...
                                      'freq', options.resample.freq, ...
                                      'read_all', 1);
% Process: Detect bad channels
bst_process('CallProcess', 'process_nst_detect_bad', sFileMocoResampled, [], ...
            'option_remove_negative', 1, ...
            'option_invalidate_paired_channels', 1, ...
            'option_max_sat_prop', 0.8);

% Convert to delta OD
redo_parent = redo_parent | options.dOD.redo;
sFile_dOD = nst_run_bst_proc([preproc_folder 'dOD'], redo_parent, 'process_nst_dOD', sFileMocoResampled, [], ...
                             'option_baseline_method', 1);  % mean
              
% High pass filter
redo_parent = redo_parent | options.high_pass_filter.redo;
sFile_dOD_filtered = nst_run_bst_proc([preproc_folder 'dOD | filtered'],  redo_parent, 'process_bandpass', sFile_dOD, [], ...
                                      'highpass', {0.005, ''}, ...
                                      'lowpass', {0, ''}, ...
                                      'attenuation', 'relax', ...
                                      'mirror', 0, ...
                                      'sensortypes', 'NIRS');

                                  
% Compute head model from full head model
redo_parent = redo_parent | options.head_model.redo;
nst_run_bst_proc([preproc_folder 'head model'], redo_parent, 'process_nst_sub_headmodel', ...
                 sFile_dOD_filtered, sFile_raw_full_head_model);                                  
                                  
% Project and convert to d[HbX]
redo_parent = redo_parent | options.projection.redo;
sFilesHbProj = nst_run_bst_proc({[preproc_folder 'dHbO_cortex'], [preproc_folder 'dHbR_cortex']},  redo_parent, ... 
                                'process_nst_cortical_projection_mne', sFile_dOD_filtered, []);

end

function [sFiles_GLM, redone_any_contrast] = glm_1st_level(sFiles, stim_events, contrasts, options, force_redo)

if nargin < 5
    force_redo = 1;
end
glm_folder = sprintf('GLM_%s/', options.PIPELINE_TAG);

[SubjectName, preprocs_folder] = bst_fileparts(bst_fileparts(sFiles{1}), 1);
sSubject = bst_get('Subject', SubjectName);

redone_any_contrast = 0; % Track if any contrast for any file had to be recomputed
sFiles_GLM = cell(1, length(sFiles));
for ifile=1:length(sFiles)
    data_cmt = load(file_fullpath(sFiles{ifile}), 'Comment');
    
    comment_glm_prefix{ifile} = [glm_folder 'GLM ' data_cmt.Comment];
    % Process: GLM - design and fit
    redo_parent = force_redo | options.GLM.redo;
    sFiles_GLM{ifile} = nst_run_bst_proc([comment_glm_prefix{ifile} ' - fit'], redo_parent, ...
                                         'process_nst_compute_glm', sFiles{ifile}, [], ...
                                         'stim_events',    strjoin(stim_events, ', '), ...
                                         'hrf_model',      1, ...  % CANONICAL
                                         'trend',          1, ...
                                         'fitting',        1, ...  % OLS
                                         'save_residuals', 0, ...
                                         'save_betas',     0);
end

for ifile=1:length(sFiles_GLM)

    for icon=1:length(contrasts)

        % Process: GLM - intra subject contrast
        redo = redo_parent | options.contrast.redo;
        %TODO: store contrast in main GLM result
        sFiles_GLM_ttest = nst_run_bst_proc([comment_glm_prefix{ifile} ' - con_t-+ ' contrasts(icon).label], redo, ...
                                            'process_nst_compute_ttest', sFiles_GLM{ifile}, [], ...
                                            'Contrast', contrasts(icon).vector, ...
                                            'tail',     'two');  
        % Plots
        data_tag = get_bst_file_tag(sFiles{ifile});
        %TODO: use current pval
        pval = 0.05;
        fig_bfn = sprintf('%s_%s_tmap_%1.1e_%s.png', SubjectName, data_tag, pval, contrasts(icon).label);
        fig_fn = protect_fn_str(fullfile(options.fig_dir, fig_bfn ));
        if options.make_figs && (redo || options.contrast_plot.redo || ~exist(fig_fn, 'file'))
            % Apply correct thresholding
            hFigSurfData = view_surface_data(sSubject.Surface(sSubject.iCortex).FileName, ...
                                             sFiles_GLM_ttest, 'NIRS', 'NewFigure');
            % TODO: set stat thresholding
            bst_figures('SetBackgroundColor', hFigSurfData, [1 1 1]);
            bst_colormaps('SetDisplayColorbar', 'stat2', 0);
            view([89 -24]); %TODO: expose as option
            zoom(hFigSurfData, 1.3);

            export_fig(fig_fn, '-transparent', '-r%d', options.fig_dpi);
            close(hFigSurfData);
        end
        
        redone_any_contrast = redone_any_contrast | redo;
    end
end


end

function file_raw_fm = get_full_head_model(sfile_raw, options)

fm_subject = ['full_head_model_' options.PIPELINE_TAG];
file_raw_fm = nst_get_bst_func_files(fm_subject, 'origin', 'Raw');

if isempty(file_raw_fm)
    warning('Subject for full head model not found.');
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
end

% Compute head model for all pairs
nst_run_bst_proc('head model [all pairs]', options.head_model.redo, ...
                 'process_nst_import_head_model', file_raw_fm, [], ...
                 'use_closest_wl', 1, 'use_all_pairs', 1);
end

function options = setup(protocol_name)

sProtocol = bst_get('Protocol', protocol_name); 
if isempty(sProtocol)
    error('Protocol %s not found', protocol_name);
end
gui_brainstorm('SetCurrentProtocol', sProtocol);

options.sci.redo = 0;

options.head_model.redo = 0;

options.moco.redo = 0;
options.moco.export_dir = create_dir(fullfile('.', 'moco_marking'));

options.resample.redo = 0;
options.resample.freq = 5; % Hz

options.dOD.redo = 0;

options.high_pass_filter.redo = 0;
options.high_pass_filter.low_cutoff = 0.01; %Hz
% TODO:
% options.high_pass_filter.high_cutoff = 0.01; %Hz

options.tag_bad_chans.redo = 0;
% TODO:
% options.tag_bad_chans.export_dir = handle_dir(fullfile('.', 'bad_chans_marking'));

options.projection.redo = 0;

options.GLM.redo = 0;

options.contrast.redo = 0;
options.contrast_plot.redo = 0;

options.make_figs = 1;
options.fig_dir = create_dir(fullfile('.', 'bad_chans_marking'));
options.fig_dpi = 100;

% Oblique view from the top
% options.plot_3d_view_az = ;
% options.plot_3d_view_el = ;

if ~function_exists('export_fig')
    disp('"export_fig" not found. Can be installed from "https://github.com/altmany/export_fig"');
    return
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
    warning('Data folder should not be part of nirstorm source folders');
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