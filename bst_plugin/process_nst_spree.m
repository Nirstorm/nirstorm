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
    
    %TODO: add option flag to rescale response to signal fit
    
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
    
   
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


