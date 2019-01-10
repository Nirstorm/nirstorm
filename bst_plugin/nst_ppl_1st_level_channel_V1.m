function varargout = nst_ppl_1st_level_channel_V1(arg1, arg2, glm_contrasts, ...
                                                  options, force_redo)
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
options.PIPELINE_TAG = '_nspc_V1';

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

if nargin < 5
    force_redo = 0;
end


if ~isempty(file_raw)
    [sFile_preprocessed, redone_preprocs, preproc_folder_raw, preproc_folder_hb] = preprocs(file_raw, options, force_redo);
else
    error('Given raw data is empty');
end

[sFiles_GLM, redone_any_contrast, glm_folder] = glm_1st_level(sFile_preprocessed, glm_events_order, glm_contrasts, options, redone_preprocs | force_redo);

if options.clean_preprocessings
    full_preproc_folder = fileparts(sFiles_preprocessed{1});
    [sStudy, iStudy] = bst_get('StudyWithCondition', full_preproc_folder);
    db_delete_studies(iStudy);
end

if nargout >= 1
    varargout{1} = sFiles_GLM;
end

if nargout >= 2
    varargout{2} = redone_any_contrast;
end

if nargout >= 3 
    varargout{3} = glm_folder;
end

if nargout >= 4 && ~options.clean_preprocessings
    varargout{4} = preproc_folder_raw;
end

if nargout >= 5 && ~options.clean_preprocessings
    varargout{5} = preproc_folder_hb;
end

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

function [sFiles_GLM, redone_any_contrast, glm_folder] = glm_1st_level(sFile, stim_events, contrasts, options, force_redo)

if nargin < 5
    force_redo = 1;
end
glm_folder = sprintf('GLM_%s/', options.PIPELINE_TAG);

[SubjectName, preprocs_folder] = bst_fileparts(bst_fileparts(sFile), 1); %#ok<ASGLU>
sSubject = bst_get('Subject', SubjectName);

redone_any_contrast = 0; % Track if any contrast for any file had to be recomputed

redo_parent = force_redo | options.GLM.redo;
    
% Process: GLM - design and fit

comment_glm_prefix = [glm_folder 'GLM '];
output_names = {[comment_glm_prefix ' | fitted model']};
if options.GLM.save_fit
    output_names{end+1} = [comment_glm_prefix ' | signal fit'];
end
[sFiles_GLM, redone_fit] = nst_run_bst_proc(output_names, redo_parent, ...
                                            'process_nst_compute_glm', sFile, [], ...
                                            'stim_events',    strjoin(stim_events, ', '), ...
                                            'hrf_model',      1, ...  % CANONICAL
                                            'trend',          1, ...
                                            'fitting',        1, ...  % OLS - precoloring
                                            'hpf_low_cutoff', options.high_pass_filter.low_cutoff, ...
                                            'trim_start', options.GLM.trim_start, ...
                                            'save_residuals', 0, ...
                                            'save_betas',     0, ...
                                            'save_fit',  options.GLM.save_fit);

if 1 %HACK
    chans_oi = {'S1D1HbO', 'S1D1HbR', 'S1D2HbO', 'S1D2HbR'};
    cond_colors = [[1 125/255 125/255] ; [125/255 125/255 1]; [125/255 1 125/255]; [.5 1 0]; [0 1 .5]];

    in_data_mat = in_bst_data(sFile);
    ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', sFile));
    glm_fit_mat = in_bst_data(sFiles_GLM{2}); % GLM signal fit
    
    trim_start_sample = round(options.GLM.trim_start / diff(in_data_mat.Time(1:2)));
    trimmed_smpl = (trim_start_sample+1):length(in_data_mat.Time);
    
    in_data_mat.F = in_data_mat.F * 1e6; %convert to mumol
    
    for ichan_oi=1:length(chans_oi)
        ichan = strcmp({ChannelMat.Channel.Name}, chans_oi{ichan_oi});
        hfig = figure(); hold on;
        ymin = min(in_data_mat.F(ichan, trimmed_smpl));
        ymax = max(in_data_mat.F(ichan, trimmed_smpl));
        plot_paradigm(in_data_mat.Events, ymin, ymax, cond_colors);
        plot(in_data_mat.Time(trimmed_smpl), in_data_mat.F(ichan, trimmed_smpl), 'b', 'LineWidth', 2);
        plot(in_data_mat.Time(trimmed_smpl), glm_fit_mat.F(ichan, trimmed_smpl), 'r', 'LineWidth', 2);
        saveas(hfig, fullfile(options.fig_dir, sprintf('%s_glm_fit_%s.png', SubjectName, chans_oi{ichan_oi})));
        close(hfig);
    end
end
                                        
for icon=1:length(contrasts)
    % Process: GLM - intra subject contrast
    redo = redone_fit | options.contrast.redo;
    %TODO: store contrast in main GLM result
    [sFiles_GLM_ttest, redone_con] = nst_run_bst_proc([comment_glm_prefix ' | con_t-+ ' contrasts(icon).label], redo, ...
                                                      'process_nst_compute_ttest', sFiles_GLM{1}, [], ...
                                                      'Contrast', contrasts(icon).vector, ...
                                                      'tail',     'two');  
    % Plots
    data_tag = get_bst_file_tag(sFiles_GLM{1});
    
    %TODO: use current pval
    pval = 0.05;

    if options.make_figs && (redo || options.contrast_plot.redo || ~exist(fig_fn, 'file'))
        % Apply correct thresholding
        hFigTopoData = view_topography(sFiles_GLM_ttest, 'NIRS', '3DOptodes');
        
        % TODO: set stat thresholding

        hb_types = {'HbO', 'HbR'}; 
        for ihb=1:length(hb_types)
            panel_montage('SetCurrentMontage', hFigTopoData, [hb_types{ihb} '[tmp]']);
            bst_figures('SetBackgroundColor', hFigTopoData, [1 1 1]);
            bst_colormaps('SetDisplayColorbar', 'stat2', 1);
            view([89 24]); %TODO: expose as option
            %zoom(hFigTopoData, 1.3);
            fig_bfn = sprintf('%s_topo_%s_tmap_%1.1e_%s.png', SubjectName, hb_types{ihb}, pval, contrasts(icon).label);
            fig_fn = protect_fn_str(fullfile(options.fig_dir, fig_bfn ));
            export_fig(fig_fn, '-transparent', sprintf('-r%d', options.fig_dpi));
        end

        close(hFigTopoData);
    end

    redone_any_contrast = redone_any_contrast | redone_con;
end
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

options.resample.redo = 0;
options.resample.freq = 5; % Hz

options.MBLL.redo = 0;
options.MBLL.age = 60;
options.MBLL.pvf = 50;
options.MBLL.baseline_method = 2; % 1:mean, 2:median

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