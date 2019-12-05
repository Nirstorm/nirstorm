function varargout = process_nst_spree( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Thomas Vincent (2016-2019)
%
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Spline-regularized response estimation';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - Wip';
    sProcess.Index       = 1409;
    sProcess.Description = '';
    

    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw', 'results'};
    sProcess.OutputTypes = {'data', 'data', 'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.stim_events.Comment = 'Stimulation events: ';
    sProcess.options.stim_events.Type    = 'text';
    sProcess.options.stim_events.Value   = '';
        
    sProcess.options.trim_start.Comment = 'Ignore starting signal: ';
    sProcess.options.trim_start.Type    = 'value';
    sProcess.options.trim_start.Value   = {0, 'sec', 2};
    
    sProcess.options.nb_iterations.Comment = 'Nb iterations';
    sProcess.options.nb_iterations.Type    = 'value';
    sProcess.options.nb_iterations.Value   = {2000, '', 0};
    
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    sProcess.options.output_cmt.Comment = '<B>Extra outputs</B>:';
    sProcess.options.output_cmt.Type    = 'label';
        
    sProcess.options.save_betas.Comment = 'Effect maps';
    sProcess.options.save_betas.Type    = 'checkbox';
    sProcess.options.save_betas.Value   =  0;
    
    sProcess.options.save_fit.Comment = 'Fit';
    sProcess.options.save_fit.Type    = 'checkbox';
    sProcess.options.save_fit.Value   =  0;
    
    sProcess.options.save_full_fitted_model.Comment = 'Save all';
    sProcess.options.save_full_fitted_model.Type    = 'checkbox';
    sProcess.options.save_full_fitted_model.Hidden  = 1;
    sProcess.options.save_full_fitted_model.Value   =  0;
    
    sProcess.options.save_fit.Comment = 'Fit';
    sProcess.options.save_fit.Type    = 'checkbox';
    sProcess.options.save_fit.Value   =  0;
    
    sProcess.options.output_fig_dir.Comment = 'Figure output directory: ';
    sProcess.options.output_fig_dir.Type    = 'text';
    sProcess.options.output_fig_dir.Hidden  = 1;
    sProcess.options.output_fig_dir.Value   = '';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
    
    save_fit = sProcess.options.save_fit.Value;

       
    %% Select events
    % TODO: utests with all events found, no event selected, no available
    %       event in data, one event missing
    if isempty(sProcess.options.stim_events.Value)
         bst_error('No event selected');
    end
    selected_event_names = cellfun(@strtrim, strsplit(sProcess.options.stim_events.Value, ','),...
                                   'UniformOutput', 0);

    % Get option values   
    nb_iterations   = sProcess.options.nb_iterations.Value{1};
    
    % export_response_figs = sProcess.options.option_do_export_response_figs.Value;
    
    DataMat = in_bst_data(sInputs(1).FileName);

    if isfield(DataMat, 'SurfaceFile')
        surface_data = 1;
        parent_data = in_bst_data(DataMat.DataFile);
        % Make sure time axis is consistent
        assert(all(parent_data.Time == DataMat.Time));
        if isempty(DataMat.Events) && isfield(parent_data, 'Events')
            DataMat.Events = parent_data.Events;
        end
        if isempty(DataMat.F) && ~isempty(DataMat.ImageGridAmp) && size(DataMat.ImageGridAmp, 2)==length(DataMat.Time)
            Y = DataMat.ImageGridAmp';
            if issparse(Y)
                Y = full(Y);
            end
        else
            bst_error('Cannot get signals from surface data');
        end
    else
        surface_data = 0;
        % Get signals of NIRS channels only:
        ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
        [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');
        Y = DataMat.F(nirs_ichans,:)';
    end
    
    % Convert to micromol
    if strcmp(DataMat.DisplayUnits, 'mol.l-1')
        Y = Y * 1e6;
        DataMat.DisplayUnits = 'mumol.l-1';
    elseif strcmp(DataMat.DisplayUnits, 'mmol.l-1')
        Y = Y * 1e3;
        DataMat.DisplayUnits = 'mumol.l-1';
    elseif ~strcmp(DataMat.DisplayUnits, 'mumol.l-1')
        if ~isempty(DataMat.DisplayUnits)
            warning('Cannot interpret data unit: %s.', DataMat.DisplayUnits);
        else
            warning('Unspecified data unit.');
        end
    end
    
    all_event_names = {DataMat.Events.label};
    events_found = ismember(selected_event_names, all_event_names);
    if ~all(events_found)
        bst_error(sprintf('Event names "%s" not found (available events: "%s")', ...
                          strjoin(selected_event_names(~events_found), ', '), ...
                          strjoin(all_event_names, ',')));
        return;
    end
    ievents = cellfun(@(l) find(strcmp(l,all_event_names)), selected_event_names);
    
    % Check that all selected events are extended
    % TODO: also handle single events for event-related design
    isExtended = false(1,length(ievents));
    for ievt=1:length(ievents)
        isExtended(ievt) = (size(DataMat.Events(ievents(ievt)).times, 1) == 2);
    end
    if ~all(isExtended)
         bst_error(sprintf('Simple events not supported: %s ', ...
                           strjoin(selected_event_names(ievents(~isExtended)), ', ')));
         return;
    end

    response_duration = 30; % seconds
    [responses, responses_std, response_time_axis, summary_stats, fmodel] = ...
        Compute(Y, DataMat.Time, DataMat.Events(ievents), response_duration, nb_iterations); 
        
%     if export_response_figs
%         fig_prefix = [sInputs.SubjectName '_' sInputs.Condition '_'];
%         plot_responses(responses, responses_std, response_time_axis, ...
%                        residuals, summary_stats, ChanneMat, fig_dir.Value{1}, fig_prefix);
%     end
    
    output_prefix = [sInputs(1).Comment ' | Spree '];
    extra_output.DisplayUnits = DataMat.DisplayUnits; %TODO: check scaling
    extra_output.spree_summary_stats = summary_stats;
    if sProcess.options.save_full_fitted_model.Value
        extra_output.spree_fmodel = fmodel;
    end
        
    output_comment = [output_prefix ' result'];
    if surface_data
        [sStudy, ResultFile] = nst_bst_add_surf_data(responses', response_time_axis, [], 'surf_spree_res', output_comment, ...
                                                     [], sStudy, 'Spree', DataMat.SurfaceFile, 0, extra_output);
        OutputFile = ResultFile;
    else
    
        nb_channels = length(ChannelMat.Channel);
        nb_response_samples = length(response_time_axis);

        responses_full = nan(nb_response_samples, nb_channels);
        responses_full(:, nirs_ichans) = responses;

        responses_std_full = nan(nb_response_samples, nb_channels);
        responses_std_full(:, nirs_ichans) = responses_std;

        sStudy = bst_get('Study', sInputs(1).iStudy); 

        % Save time-series data
        %sDataOut = db_template('data');
        sDataOut.F            = responses_full';
        sDataOut.Comment      = output_comment;
        sDataOut.ChannelFlag  = DataMat.ChannelFlag; %ones(size(responses, 2), 1);
        sDataOut.Time         = response_time_axis';
        sDataOut.DataType     = 'recordings'; 
        sDataOut.nAvg         = 1;
        sDataOut.Std = responses_std_full';
        sDataOut.ColormapType = [];
        sDataOut.Events = [];

        % Add extra fields
        extra_fields = fieldnames(extra_output);
        for ifield = 1:length(extra_fields)
            sDataOut.(extra_fields{ifield}) = extra_output.(extra_fields{ifield});
        end
        
        % Generate a new file name in the same folder
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),...
                                 'data_spree');
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs(1).iStudy, OutputFile, sDataOut);
        
        
        
        ChannelMat
    end
    
    if save_fit
        fit = fmodel.signal_fit;
        output_tag = sprintf('ir%d_spree_signal_fit', sInputs(1).iItem);
        output_comment = [output_prefix '- signal fit'];
        if surface_data
            [sStudy, ResultFile] = nst_bst_add_surf_data(fit', DataMat.Time, [], output_tag, output_comment, ...
                                                         [], sStudy, 'Spree', DataMat.SurfaceFile);
        else
            data_out = zeros(size(DataMat.F));
            data_out(nirs_ichans, :) = fit';
            Out_DataMat = db_template('data');
            Out_DataMat.F           = data_out;
            Out_DataMat.Comment     = output_comment;
            Out_DataMat.DataType     = 'recordings';
            Out_DataMat.Time        =  DataMat.Time;
            Out_DataMat.Events      =  DataMat.Events;
            Out_DataMat.ChannelFlag =  DataMat.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
            Out_DataMat.DisplayUnits = DataMat.DisplayUnits;
            Out_DataMat.nAvg         = 1;
            
            Out_DataMat = bst_history('add', Out_DataMat, DataMat.History, '');
            Out_DataMat = bst_history('add', Out_DataMat, 'Spree', FormatComment(sProcess));
            fit_fn = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'data_spree_fit');
            Out_DataMat.FileName = file_short(fit_fn);
            bst_save(fit_fn, Out_DataMat, 'v7');
            db_add_data(sInputs(1).iStudy, fit_fn, Out_DataMat);
        end        
    end
    
    fig_dir = sProcess.options.output_fig_dir.Value;
    if ~isempty(fig_dir)
        if ~exist(fig_dir, 'dir')
            mkdir(fig_dir);
        end
        fig_prefix = '';
        close_figs = 1; 
        plot_deconv_results(fmodel, {ChannelMat.Channel(nirs_ichans).Name}, fig_dir, fig_prefix, close_figs)
    end
    
    OutputFile = {OutputFile};
end

function [responses, responses_pstd, response_time_axis, summary_stats, fitted_model] = ...
    Compute(nirs_sig, time, events, response_duration, nb_iterations)
%% Apply 
% Args
%    - nirs_sig: matrix of double, size: nb_samples x nb_channels
%        Measured nirs Hb signal in micromol.L-1.cm-1.
%    - time: array of double, size: nb_samples x 1
%        Time axis of the input nirs signal.
%    - events:
%    - response_duration: float
%      Response duration in second.
%    [- nb_iterations ]: positive double, default is 50
% 
% Output: 
%   - responses: matrix of double, size: nb_samples x nb_channels
%     Responses in delta [HbO/R/T] in micromol.L-1.cm-1.
%     These are "block" responses, ie the convolution of the average stimulus with the HRFs. 
%     This is to have a more meaningful scale that the HRF which associated 
%     with a Dirac stimulus event.
%   - responses_std: matrix of double, size: nb_samples x nb_channels
%     Response standard deviations in delta [HbO/R/T] in micromol.L-1.cm-1 
%   - response_time_axis: 1d array of double
%   - summary_stats: struct with fields:
%         * activating_positions
%         * + PPM stats
%         * time_to_peak
%         * time_to_peak_std
%     

if nargin < 5
    nb_iterations = 50;
end

nirs.data = nirs_sig;
nirs.time = time;

if length(events) > 1
    warning('Only one condition is handled. Taking only 1st one');
end

nirs.events = events(1);
% for ievt=1:length(events)
%     nirs.paradigm.conditions{ievt} = events(ievt).label;
%     nirs.paradigm.onsets{ievt} = events(ievt).times(1,:);
%     if size(events(ievt).times, 1) == 2 % extended events
%         nirs.paradigm.durations{ievt} = diff(events(ievt).times, 1);
%     else % simple events -> assume brief stimuli
%         nirs.paradigm.durations{ievt} = diff(time(1:2));
%     end
% end

%% Setup model
model = spree_splhblf_init_default(); 
model.options.trend_type = 'cosine';
model.options.trend_bands(1).name = 'baseline';
model.options.trend_bands(1).bounds = [0 1/response_duration]; %[0 0.01]; % .[0 0.0225]
model.options.trend_bands(2).name = 'respiration';
model.options.trend_bands(2).bounds = [1/5 1/3];
model.options.trend_bands(3).name = 'heart';
model.options.trend_bands(3).bounds = [0.8 1.7];
model.options.npb_mirror_trend_fix = 150; %nb  samples to mirror at boundaries

model.options.knots_type = 'regular';
model.options.response_duration = 30;
model.options.regular_knots_dt = 5;
model.options.bandpass_residuals = [];
model = spree_set_estimation(model, {'f1', 'noise_var', 'f1_var', 'response', 'trend_coeffs', 'trend_var'}); % model_set_estimation
model.max_iterations = nb_iterations;
model.var_hist_pace = 1; %TODO: debug when pace > 1
model.obs_hist_pace = 1;
model = spree_splhblf_init_from_data(model, nirs);
model.burnin_period = round(model.max_iterations/3);

%% Fit model
fitted_model = spree_splhblf_fit(model);

%% Outputs
responses = fitted_model.observables.response_block_pm;
responses_pstd = fitted_model.observables.response_block_pstd;
response_time_axis = fitted_model.constants.response_block_time_axis;
summary_stats = fitted_model.summary;
end


function [fdata, fchannel_def, data_other, channel_def_other] = ...
    filter_data_by_channel_type(data, channel_def, channel_types)

%    - channel_types: str or cell array of str
%        Channel types to keep.
%        If [] is given then all channels are kept.

if nargin < 3 || isempty(channel_types)
    channel_types = unique({channel_def.Channel.Type});
else
    if isstr(channel_types)
        channel_types = {channel_types};
    end
end

kept_ichans = ismember({channel_def.Channel.Type}, channel_types);
fchannel_def = channel_def;
fchannel_def.Channel = channel_def.Channel(kept_ichans);
fdata = data(:, kept_ichans);

other_ichans = ~kept_ichans;
channel_def_other = channel_def;
channel_def_other.Channel = channel_def.Channel(other_ichans);
data_other = data(:, other_ichans);
end


function [pair_names, pair_indexes] = get_pairs(channel_def)
%% group paired channels -> explode channel data according to pairs
% Args
%    - channel_def: struct
%        Defintion of channels as given by brainstorm
%        Used fields: Nirs.Hb, Channel
%
% ASSUME: data contain only Hb-related channels (no AUX etc.)
%
% TOCHECK WARNING: uses containers. Map which is available with matlab > v2008
%
%  Outputs: 
%     - pair_names: cell array of str, size: nb_pairs
%         Pair names, format: SXDX
%     - pair_indexes: matrix of double, size: nb_pairs x nb_hbs
%         Input channel indexes grouped by pairs
%

nb_hbs = length(channel_def.Nirs.Hb);

pair_to_chans = containers.Map();
for ichan=1:length(channel_def.Channel)
    chan_name = channel_def.Channel(ichan).Name;
    ihb = strfind(chan_name, 'Hb'); 
    pair_name = chan_name(1:ihb-1);
    hb = chan_name(ihb:end);
    if pair_to_chans.isKey(pair_name)
        hba = pair_to_chans(pair_name);
    else
        hba = zeros(1, nb_hbs);
    end
    hba(strcmp(hb, channel_def.Nirs.Hb)) = ichan;
    pair_to_chans(pair_name) = hba;
end
nb_pairs = length(channel_def.Channel) / nb_hbs;
pair_names = pair_to_chans.keys;
pair_indexes = zeros(nb_pairs, nb_hbs);
for ipair=1:nb_pairs
    p_indexes = pair_to_chans(pair_names{ipair});
    pair_indexes(ipair, :) = p_indexes;
end

end


function plot_deconv_results(model, chan_labels, output_dir, fig_prefix, close_figs)

save_figs = 1;

peak_types = {'peak', 'undershoot'};
extremas = {'max', 'min'};
peak_flag = [2 1];
tt_ic_colors = {[.4, 0, 0], [0, 0, .4]};
curve_colors = {[1, .2, 0],  [0, .2, 1]};
pbar_colors = {[.6 0 0], [0 0 .6]};
new_ylim_top = nan;
new_ylim_bot = nan;
gmax = nan;
gmin = nan;

% hrf_ylims = [-0.8 0.6];
% cmro2_cbf_ylims =  [0.6 2];


%% save summary stats + effects to csv file
ic = 1;
chan_names = {};

col_names = {'peak_size', 'peak_time', 'peak_IC_low', 'peak_IC_up', 'peak_post_prob', 'dip_size', 'dip_time', 'dip_IC_low', 'dip_IC_up', 'dip_post_prob', 'extremum_size', 'extremum_pval'};
for ichan=1:model.constants.n_channels

    chan_label = chan_labels{ichan};
    
    chan_names{ic} = chan_label;
        
    peak_size(ic) = model.summary.a_max_pm(ichan);
    peak_time(ic) = model.summary.t_a_max_pm(ichan);
    peak_IC_low(ic) = model.summary.t_a_max_IC(ichan, 1);
    peak_IC_up(ic) = model.summary.t_a_max_IC(ichan, 2);
    peak_proba(ic) = model.summary.prob_a_max(ichan);
        
    dip_size(ic) = model.summary.a_min_pm(ichan);
    dip_time(ic) = model.summary.t_a_min_pm(ichan);
    dip_IC_low(ic) = model.summary.t_a_min_IC(ichan, 1);
    dip_IC_up(ic) = model.summary.t_a_min_IC(ichan, 2);
    dip_proba(ic) = model.summary.prob_a_min(ichan);
    
    if abs(peak_size(ic)) > abs(dip_size(ic))
        extremum_size(ic) = peak_size(ic);
        extremum_post_prob(ic) = 1-peak_proba(ic);
    else
        extremum_size(ic) = dip_size(ic);
        extremum_post_prob(ic) = 1-dip_proba(ic);
    end
    %         for ipeak_type=1:length(peak_types)
    %             plab = peak_types{ipeak_type};
    %             tt_lower = model.summary.(['tt' plab(1) '_lower'])(ichan);
    %             tt_upper = model.summary.(['tt' plab(1) '_upper'])(ichan);
    %             significance = prob_to_star(1-model.summary.(['pk_' plab '_prob_tau_hat'])(ichan));
    %             fprintf(freport_tt, [hb_type '-' plab ', %1.2f, %1.2f, %s\n'], tt_lower, tt_upper, significance);
    %         end
    ic = ic + 1;
end
summary_table = table(peak_size', peak_time', peak_IC_low', peak_IC_up', peak_proba', ...
    dip_size', dip_time', dip_IC_low', dip_IC_up', dip_proba', ...
    extremum_size', extremum_post_prob', ...
    'RowNames', chan_names, 'VariableNames', col_names);
summary_fn = fullfile(output_dir, [fig_prefix 'summary_stats_effects.csv']);
writetable(summary_table, summary_fn, 'Delimiter',';', 'WriteRowNames', 1);

for ichan=1:model.constants.n_channels
    
    chan_label = chan_labels{ichan};
    
    %% Plot responses with peaks & error bars
    flabel = [fig_prefix 'response_PM_peak_results_' chan_label];
    h = fig_from_label(flabel); hold on;
    peak_type = 'extremum';

    if ~isempty(strfind('HbO', chan_label))
        ihb = 1;
    else
        ihb = 2;
    end
    plot_response_density(model, ichan, curve_colors{ihb}, peak_type);
    plot_response_profile(model, ichan, curve_colors{ihb});
    
    font_opts.fontsize = 40;
    font_opts.rheight = 1/10;
    
    plot_ic_annotations(model, ichan, curve_colors{ihb}, peak_type, font_opts);
    
    plot_ttp_annotations(model, ichan, curve_colors{ihb}, peak_type, font_opts);
        
    % plot noise 2.5 percentiles
    %     color = [0.9 0.9 0.9];
    %     xlims = xlim();
    %     ylims = ylim();
    %     noise_1st_pctl = 1;
    %     rectangle('Position', [xlims(1) ylims(1) diff(xlims) diff(ylims)], ...
    %               'FaceColor', color, 'Linestyle', 'none');
    %
    %     hstack = get(gca, 'Children');
    %     set(gca, 'Children', [hstack(2:end) ; hstack(1)]);resid_colors
    if 0
        % plot std dev of bandpassed residuals
        resid_colors = {[255 216 207] / 255,  [207 216 255] / 255};
        for ihb=1:length(hb_types)
            model = models{ihb};
            resid_color = resid_colors{ihb};
            res = model.constants.y(:, ichan) - model.stim_induced_fit(:, ichan);
            
            sfit = model.constants.XP * model.observables.f1_pm(:, ichan);
            response_pm = model.observables.response_pm(:, ichan);
            
            if 0
                %% detrend & filter data
                bp_res = detrend(res);
                opts.ftype = 'bandpass';
                opts.low_cutoff = model.options.bandpass_residuals(1);
                opts.high_cutoff = model.options.bandpass_residuals(2);
                opts.order = 3;
                opts.flag_keepMean = 1;
                bp_res = mfip_IIR_butter(bp_res, model.constants.sampling_freq, opts);
                
                scale_factor = (max(response_pm) - min(response_pm)) / (max(sfit) - min(sfit));
                bp_res_rescaled = bp_res * scale_factor;
                scale_factors{ihb} = scale_factor;
                pctl = 5; % percentage
                npctls = round(1/pctl * 100);
                qtls = quantile(bp_res_rescaled, npctls);
                % qtls = [min(bp_res_rescaled) max(bp_res_rescaled)];
                xlims = xlim();
                rectangle('Position', [xlims(1) qtls(1) diff(xlims) qtls(end)-qtls(1)], ...
                    'FaceColor', resid_color, 'Linestyle', 'none');
            else
                scale_factors{ihb} = 1;
            end
            if 0
                figure(); hold on;
                plot(model.constants.y_time_axis, model.constants.y(:, ichan), 'k')
                plot(model.constants.y_time_axis, model.signal_fit(:, ichan), 'r')
                plot(model.constants.y_time_axis, res, 'b');
                plot(model.constants.y_time_axis, sfit, 'g');
            end
        end
    end
    %hstack = get(gca, 'Children');
    %set(gca, 'Children', [hstack(3:end) ; hstack(1:2)]);
    set(gca, 'YColor', [0. 0. 0.]);
    xlim([model.constants.response_time_axis(1), model.constants.response_time_axis(end)]);
    %ylabel('$\Delta$ [Hb]');
    %xlabel('Time in sec.');
    %fill_fig(h);
    
    %plot_2d_contour_hist(model.vars_history.response_time_to_peak', model.vars_history.response_peak_value')
    if save_figs
        save_fig(flabel, output_dir);
    end
    if close_figs
        close(h);
    end
    
    if 0
        %% plot of fit vs input signal
        for ihb = 1:length(hb_types) % channel of interest
            hb_type = hb_types{ihb};
            chan_oi = chans_oi{ihb}(ichan);
            model = models{ihb};
            y_channels = model.constants.y;
            evoked_samples_oi = 1:length(model.constants.y_time_axis); % nb of samples of interest for plots
            time_axis = model.constants.y_time_axis;
            sampled_paradigm = model.constants.X(:,1);
            baseline = model.constants.T * model.observables.trend_coeffs_pm(:, ichan);
            stim_induced_fit = model.constants.XP * ...
                model.observables.f1_pm(:, ichan) + baseline;
            
            flabel = [fig_prefix '_deconv_signal_fit_' chan_label '_' hb_type];
            h = fig_from_label(flabel); hold on;
            
            plots = [];
            p = plot(time_axis(evoked_samples_oi), y_channels(evoked_samples_oi, ichan), ...
                'DisplayName', sprintf('signal in channel %d', ichan), ...
                'Color', [0 0.45 0.74], 'Linewidth', 1.5);
            plots = [plots p];
            
            p = plot(time_axis(evoked_samples_oi), baseline(evoked_samples_oi), ...
                'DisplayName', sprintf('baseline in channel %d', ichan), ...
                'Color', 'k', 'Linewidth', 1.5);
            plots = [plots p];
            
            
            p = plot(time_axis(evoked_samples_oi), stim_induced_fit(evoked_samples_oi), ...
                'DisplayName', sprintf('stim-induced fit in channel %d', ichan), ...
                'Color', 'r', 'Linewidth', 3);
            plots = [plots p];
            
            %         p = plot(time_axis(evoked_samples_oi), sampled_paradigm(evoked_samples_oi) * max(max(y_channels(evoked_samples_oi, coi))) * .2,...
            %             'Color', [1 0.5 0], 'DisplayName', sprintf('stimulus'));
            %         plots = [plots p];
            % legend(plots);
            xlim([time_axis(evoked_samples_oi(1)) time_axis(evoked_samples_oi(end))]);
            nst_plot_paradigm(events, nan);
            set(get(h, 'CurrentAxes'), 'FontSize', 16);
            set(gca,'LooseInset',get(gca,'TightInset'));
            ylabel(['$\Delta$ [' hb_types{ihb} ']']);
            %xlabel('Time (sec.)');
            if save_figs
                save_fig(flabel, output_dir);
            end
            if close_figs
                close(h);
            end
        end
    end
end


end

function plot_response_density(model, ichan, curve_color, peak_def)
color_breaks = [[100 65 65]; ...
                   [255 255 220]] / 255;          
% color_breaks = [curve_color * 0.3; [255 255 220]/255];
cmap_density = create_linear_segmented_cmap(color_breaks);
warm_idx = model.burnin_period:model.max_iterations;

if strcmp(peak_def, 'max') || ...
        (strcmp(peak_def, 'extremum') && abs(model.summary.a_max_pm(ichan)) > abs(model.summary.a_min_pm(ichan)))
    peak_ttp_samples = model.vars_history.response_time_to_peak(warm_idx, ichan);
    peak_samples = model.vars_history.response_max(warm_idx, ichan);
else
    peak_ttp_samples = model.vars_history.response_time_to_undershoot(warm_idx, ichan);
    peak_samples = model.vars_history.response_min(warm_idx, ichan);
end

alpha = 0.0001;
percentiles = [alpha/2 1-alpha/2] .* 100;
peak_ic = prctile(peak_samples, percentiles);

if var(peak_ttp_samples) > eps && var(peak_samples) > eps
    resp_max = max(max(model.vars_history.response(warm_idx, :, ichan))) * 1.01;
    resp_min = min(min(model.vars_history.response(warm_idx, :, ichan))) * 1.01;
    [f, xi] = ksdensity([peak_ttp_samples peak_samples], ...
                        'support', [model.constants.response_time_axis(1)-eps resp_min; ...
                                    model.constants.response_time_axis(end)+eps resp_max]);
    plot_contour_ksdensity(xi, f, cmap_density);

    ylims = ylim();
    ylim([min(ylims(1), peak_ic(1)), max(ylims(2), peak_ic(2))]); 
end

end

function plot_response_profile(model, ichan, curve_color)
linewdith = 4.5;

response_pm = model.observables.response_pm;
response_pstd = model.observables.response_pstd;
t_axis = model.constants.response_time_axis;

plot(t_axis, response_pm(:, ichan), 'Color', curve_color, ...
    'LineWidth', linewdith);
errorbar(t_axis(1:8:end), response_pm(1:8:end, ichan), ...
    response_pstd(1:8:end, ichan), 'Color', curve_color, ...
    'LineWidth', linewdith * 0.3, 'LineStyle', 'none');
ymin_estim = min(response_pm(:, ichan)) - max(response_pstd(:, ichan));
ymax_estim = max(response_pm(:, ichan)) + max(response_pstd(:, ichan));
ylims = ylim();
ylim([min([ylims(1), ymin_estim]), max([ylims(2), ymax_estim])]);

xlim([t_axis(1) t_axis(end)]);

end

function plot_ic_annotations(model, ichan, curve_color, peak_def, font_opts)
global fig_options;

warm_idx = model.burnin_period:model.max_iterations;

if strcmp(peak_def, 'max') || ...
        (strcmp(peak_def, 'extremum') && abs(model.summary.a_max_pm(ichan)) > abs(model.summary.a_min_pm(ichan)))
    peak_ttp_samples = model.vars_history.response_time_to_peak(warm_idx, ichan);
    peak_samples = model.vars_history.response_max(warm_idx, ichan);
    a_pm = model.summary.a_max_pm(ichan);
    t_a_pm = model.summary.t_a_max_pm(ichan);
    t_a_IC = model.summary.t_a_max_IC(ichan, :);
    f_pm = model.summary.max_f_pm(ichan);
    peak_is_max = 1;
else
    peak_ttp_samples = model.vars_history.response_time_to_undershoot(warm_idx, ichan);
    peak_samples = model.vars_history.response_min(warm_idx, ichan);
    a_pm = model.summary.a_min_pm(ichan);
    t_a_pm = model.summary.t_a_min_pm(ichan);
    t_a_IC = model.summary.t_a_min_IC(ichan, :);
    f_pm = model.summary.min_f_pm(ichan);
    peak_is_max = 0;
end

alpha = 0.01;
percentiles = [alpha/2 1-alpha/2] .* 100;
peak_ic = prctile(peak_samples, percentiles);

ylims = ylim();
dy = diff(ylims);
t_lims = xlim();

txt_fontsize = font_opts.fontsize;
txt_rheight = font_opts.rheight;
txt_height = dy * txt_rheight;
txt_ci_width = 13; %sec

if peak_is_max
    y_txt_ci = ylims(2);
    ylim([ylims(1), ylims(2) + txt_height]);
    txt_valign = 'bottom';
else
    y_txt_ci = ylims(1);
    ylim([ylims(1) - txt_height, ylims(2)]);
    txt_valign = 'top';
end

if t_a_IC(2) <= t_lims(2) - txt_ci_width
    x_txt_ci_left = t_a_IC(2);
    x_txt_ci_right = x_txt_ci_left + txt_ci_width;
    x_txt_ci = x_txt_ci_left;
    txt_halign = 'left';
    x_txt_ci = x_txt_ci_left;
elseif t_a_IC(1) >= txt_ci_width
    x_txt_ci_right = t_a_IC(1);
    x_txt_ci_left = x_txt_ci_right - txt_ci_width;
    txt_halign = 'right';
    x_txt_ci = x_txt_ci_right;
else
    x_txt_ci = mean(t_lims);
    txt_halign = 'center';
    x_txt_ci_left = x_txt_ci - txt_ci_width/2;
    x_txt_ci_right = x_txt_ci + txt_ci_width/2;
end

linewdith = 3;
line([t_a_IC(1) x_txt_ci_left], ...
     [a_pm y_txt_ci], ...
    'Linestyle', ':', 'Color', curve_color * 0.6, 'Linewidth', linewdith * 0.75);
line([t_a_IC(2) x_txt_ci_right], ...
     [a_pm y_txt_ci], ...
    'Linestyle', ':', 'Color', curve_color * 0.6, 'Linewidth', linewdith * 0.75);

txt_ic = ['$$' sprintf('\\textrm{CI}_{\\tau,0.95}=[%s;%s]s', ...
                       pad(t_a_IC(1)), pad(t_a_IC(2))) '$$'];
                   
htxt_a = text(x_txt_ci, y_txt_ci, txt_ic, ...
            'VerticalAlignment', txt_valign, ...
            'HorizontalAlignment', txt_halign, 'fontsize', txt_fontsize * fig_options.res_factor, ...
            'BackgroundColor', 'white', ...
            'interpreter', 'latex', 'Color',  curve_color * 0.7);

end

function s = pad(n)
if n >= 10
    s = sprintf('%1.1f', n);
else
    s = sprintf('\\,%1.1f', n);
end

end

function plot_ttp_annotations(model, ichan, curve_color, peak_def, font_opts)
global fig_options;

ylims = ylim();
dy = diff(ylims);
txt_fontsize = font_opts.fontsize;
txt_rheight = font_opts.rheight;
txt_height = dy * txt_rheight;
txt_width = 7.25;
vgap = dy * 0.05;
ylim([ylims(1)-txt_height-vgap, ylims(2)]);

linewdith = 4;
y_txt = ylims(1) - vgap;


if strcmp(peak_def, 'max') || ...
        (strcmp(peak_def, 'extremum') && abs(model.summary.a_max_pm(ichan)) > abs(model.summary.a_min_pm(ichan)))
    peak = model.summary.max_f_pm(ichan);
    tpeak = model.summary.t_a_max_pm(ichan);
    peak_is_max = 1;
else
    peak = model.summary.min_f_pm(ichan);
    tpeak = model.summary.t_a_min_pm(ichan);
    peak_is_max = 0;
end
    
line([tpeak, tpeak], ...
     [y_txt, peak], ...
    'Linewidth', linewdith * 0.75, 'Color',  curve_color, 'LineStyle', '--');

x_txts = tpeak;

txt_halign = 'center';

text(x_txts, y_txt, ...
    ['$$' sprintf('\\widehat{\\tau}=%1.1f s', tpeak) '$$'], 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', txt_halign, 'fontsize', txt_fontsize * fig_options.res_factor, ...
    'Color', curve_color, 'interpreter', 'latex');

end

function cmap = create_linear_segmented_cmap(color_breaks, cmap_size, stops)

if nargin < 2
    cmap_size = 256;
end

if nargin < 3
    stops = linspace(0, 1, size(color_breaks, 1));
end

cmap_values = linspace(0, 1, cmap_size);
cmap = interp1(stops, color_breaks, cmap_values);

if 0
    figure(); imagesc(1:(min(100, cmap_size)));
    colormap(cmap);
end

end


function dminy = plot_contour_ksdensity(xi, f, cmap)

% cmap = colormap('pink');
% cmap(1,:) = [80, 50, 50] / 255;

if nargin < 3
    color_breaks = [[100 65 65]; ...
        [255 255 220]] / 255;
    cmap = create_linear_segmented_cmap(color_breaks);
end

% hh = hist(f, size(cmap,1));
% f(f < hh(end)) = 0;
x1 = xi(:,1); 
x2 = xi(:,2);
fnorm = f / max(f);
[xq,yq,z] = computeGrid(x1,x2,fnorm);
colormap(cmap);
[C, h] = contourf(xq,yq,z, 0.025:0.05:1, 'LineStyle', 'none');
dminy = min(yq(:));
end

function [xq,yq,z] = computeGrid(x1,x2,fout)
x = linspace(min(x1),max(x1));
y = linspace(min(x2),max(x2));
[xq,yq] = meshgrid(x,y);
orig_state = warning;
warning('off','all');
z = griddata(x1,x2,fout,xq,yq);
warning(orig_state);
end