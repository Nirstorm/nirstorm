function varargout = process_nst_glm_fit( varargin )
% process_compute_glm: compute the glm : find B such as Y = XB +e with X
% 
% OlS_fit use an ordinary least square algorithm to find B : B= ( X^{T}X)^{-1} X^{T} Y 
% AR-IRLS : Details about AR-IRLS algorithm can be found here :
%   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3756568/ 
% 
% 
% Further update : use more sophisticated method to fit the B
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
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
% Authors: Edouard Delaire, Thomas Vincent 2018-2019
%
% TODO: use different output types for channel-space and surface analyses
%       -> more convenient for contrast computation afterwards, to be able to map
%       input to output types
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - 1st level design and fit';
    sProcess.Category    = 'File2';
    sProcess.SubGroup    = {'NIRS', 'GLM'};
    sProcess.Index       = 1601;
    sProcess.isSeparator = 0;

    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw', 'results'};
    sProcess.OutputTypes = {'data', 'data', 'results'};
    
    sProcess.nInputs     = 2;
    sProcess.isPaired    = 0; % Should it be 1??
    sProcess.nMinFiles   = 1;
    sProcess.nOutputs    = 1;
    
    sProcess.options.label0.Comment = '<U><B>Signal Information</B></U>:';
    sProcess.options.label0.Type    = 'label';
    
    sProcess.options.filter_model.Comment = {'No filter', 'IIR filter','FIR filter', 'Filter applied on the data: '; ...
                                            'none', 'IIR_bp', 'FIR_bp',''};
    sProcess.options.filter_model.Type    = 'radio_linelabel';
    sProcess.options.filter_model.Value   = 'IIR_bp';
    sProcess.options.filter_model.Controller = struct('IIR_bp','IIR_bp','FIR_bp','FIR_bp' );
    
    sProcess.options.hpf_low_cutoff.Comment = 'High-pass filter frequency: ';
    sProcess.options.hpf_low_cutoff.Type    = 'value';
    sProcess.options.hpf_low_cutoff.Value   = {0.01, 'Hz', 2};

    sProcess.options.hpf_high_cutoff.Comment = 'High-pass filter frequency: ';
    sProcess.options.hpf_high_cutoff.Type    = 'value';
    sProcess.options.hpf_high_cutoff.Value   = {0.1, 'Hz', 2};

    sProcess.options.hpf_transition_band.Comment = 'Transition band (0=default):';
    sProcess.options.hpf_transition_band.Type    = 'value';
    sProcess.options.hpf_transition_band.Value   = {0.05, 'Hz', 2};
    sProcess.options.hpf_transition_band.Class = 'FIR_bp';

    sProcess.options.hpf_order.Comment = 'Filter order:';
    sProcess.options.hpf_order.Type    = 'value';
    sProcess.options.hpf_order.Value   = {3, '', 0};
    sProcess.options.hpf_order.Class = 'IIR_bp';


    sProcess.options.has_detrending.Comment = 'Did you used NIRS > remove slow-fluctuations ? ';
    sProcess.options.has_detrending.Type    = 'checkbox';
    sProcess.options.has_detrending.Value   =  0;
    
    sProcess.options.dct_cutoff.Comment = 'Detrend miminum Period (0= only linear detrend): ';
    sProcess.options.dct_cutoff.Type    = 'value';
    sProcess.options.dct_cutoff.Value   = {200, 's', 0};
        
    sProcess.options.label3.Comment = '<U><B>Design Matrix</B></U>:';
    sProcess.options.label3.Type    = 'label';
    
    sProcess.options.stim_events.Comment = 'Stimulation events: ';
    sProcess.options.stim_events.Type    = 'text';
    sProcess.options.stim_events.Value   = '';
    
    sProcess.options.hrf_model.Comment = 'HRF model: ';
    sProcess.options.hrf_model.Type    = 'combobox';
    sProcess.options.hrf_model.Value   = {1, fieldnames(get_hrf_types())};
    
    sProcess.options.label4.Comment = '<U><B>Nuisance Regressor</B></U>:';
    sProcess.options.label4.Type    = 'label';
    
    sProcess.options.lfO.Type    = 'radio_line';
    sProcess.options.lfO.Comment   = {'constant','constant+Linear trend','constant+Linear trend+DCT','Slow fuctuations: '};
    sProcess.options.lfO.Value=1;
    
    sProcess.options.lfo_cutoff.Comment = 'Detrend miminum Period: ';
    sProcess.options.lfo_cutoff.Type    = 'value';
    sProcess.options.lfo_cutoff.Value   = {200, 's', 0};

    % === FREQ BANDS
    %sProcess.options.use_DCT.Comment = 'Adittionals DCT';
    %sProcess.options.use_DCT.Type    = 'checkbox';
    %sProcess.options.use_DCT.Value   =  0;
 
    %sProcess.options.freqbands.Comment = '';
    %sProcess.options.freqbands.Type    = 'groupbands';
    %sProcess.options.freqbands.Value   =  getDefaultFreqBands();
    
    sProcess.options.SS_chan.Type    = 'radio_line';
    sProcess.options.SS_chan.Comment   = {'No','Based on Source-Detector distances','Based on Names','Short-separation channels: '};
    sProcess.options.SS_chan.Value=0;
    
    sProcess.options.SS_chan_distance.Comment = 'Distance threeshold: ';
    sProcess.options.SS_chan_distance.Type    = 'value';
    sProcess.options.SS_chan_distance.Value   = {1.5, 'cm', 1};
    
    sProcess.options.SS_chan_name.Comment = 'Superfical Channel [coma-separated list]';
    sProcess.options.SS_chan_name.Type    = 'text';
    sProcess.options.SS_chan_name.Value   = '';     
    
    
    sProcess.options.label2.Comment = '<U><B>Serial Correlation Preprocessing</B></U>:';
    sProcess.options.label2.Type    = 'label';
    
    sProcess.options.statistical_processing.Type    = 'radio_line';
    sProcess.options.statistical_processing.Comment   = {'Pre-coloring', 'Pre-whitenning','Method : '};
    sProcess.options.statistical_processing.Value   = 1;
        
    
    sProcess.options.output_cmt.Comment = '<U><B>Extra outputs</B></U>:';
    sProcess.options.output_cmt.Type    = 'label';
    
    sProcess.options.extra_output.Comment = 'Save extra outputs (B maps, residuals, and fit)';
    sProcess.options.extra_output.Type    = 'checkbox';
    sProcess.options.extra_output.Value   =  0;
    
    sProcess.options.save_betas.Comment = 'Beta maps';
    sProcess.options.save_betas.Type    = 'checkbox';
    sProcess.options.save_betas.Value   =  0;
    sProcess.options.save_betas.Hidden   = 1;
    
    sProcess.options.save_residuals.Comment = 'Residuals';
    sProcess.options.save_residuals.Type    = 'checkbox';
    sProcess.options.save_residuals.Value   =  0;
    sProcess.options.save_residuals.Hidden   = 1;
     
    sProcess.options.save_fit.Comment = 'Fit';
    sProcess.options.save_fit.Type    = 'checkbox';
    sProcess.options.save_fit.Value   =  0;
    sProcess.options.save_fit.Hidden   = 1;
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment; 
end

function OutputFiles = Run(sProcess, sInput, sInput_ext) %#ok<DEFNU>

    OutputFiles = {};

    basis_choice = sProcess.options.hrf_model.Value{1};

    save_residuals = sProcess.options.extra_output.Value || sProcess.options.save_residuals.Value;
    save_betas = sProcess.options.extra_output.Value ||  sProcess.options.save_betas.Value;
    save_fit = sProcess.options.extra_output.Value || sProcess.options.save_fit.Value;

    %% Select events
    if isempty(sProcess.options.stim_events.Value)
         bst_error('No event selected');
         return;
    end
    selected_event_names = cellfun(@strtrim, strsplit(sProcess.options.stim_events.Value, ','),...
                                   'UniformOutput', 0);
    
    % Load Channels informations
    ChannelMat = in_bst_channel(sInput.ChannelFile);

    % Load recordings
    surface_data = 0;

    if strcmp(sInput.FileType, 'raw')  % Continuous data file       
        DataMat = in_bst(sInput.FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInput.FileName, 'F');
        DataMat.Events = sDataRaw.F.events;

    elseif strcmp(sInput.FileType, 'data')     % Imported data structure
        DataMat = in_bst_data(sInput.FileName);

    elseif strcmp(sInput.FileType, 'results')  % Imported data on the cortex
        surface_data = 1;
        
        DataMat = in_bst_data(sInput.FileName);
        channel_data = in_bst_data(DataMat.DataFile);

        % Make sure time axis is consistent
        assert(all(channel_data.Time == DataMat.Time));
        if isempty(DataMat.Events) && isfield(channel_data, 'Events')
            DataMat.Events = channel_data.Events;
        end
    end

    data_types = {}; % Chromophore present in the data (eg HbO, HbR, or WL860)
    if surface_data
        if ~isempty(channel_data.F) && ~isempty(DataMat.ImageGridAmp) && size(DataMat.ImageGridAmp, 2)==length(DataMat.Time)
            Y = DataMat.ImageGridAmp';
            if issparse(Y)
                Y = full(Y);
            end

            n_voxel = size(Y,2);
            mask    = find(~all(Y==0,1)); % Only keep channel in the field of view
            Y       = Y(:,mask);

        else
            bst_error('Cannot get signals from surface data');
        end
        
        hb_types = {'hbo','hbr','hbt','wl'};
        iChromohpore = find(cellfun(@(x)contains(lower(sInput.FileName),x),hb_types));
        if iChromohpore < 4
            data_types = hb_types(iChromohpore);
        else
            data_types = regexp(lower(sInput.FileName),'wl([0-9]*)','match'); %note: there is probably an easier way but it works :)
        end

    else

        % Get signals of NIRS channels only:
        [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');
        Y       = DataMat.F(nirs_ichans,:)';
        mask    = [];
        n_voxel = size(Y,2);
        data_types = unique({ChannelMat.Channel.Group});
    end
    
    % TODO: Check for dOD
    Y=nst_misc_convert_to_mumol(Y,DataMat.DisplayUnits);
    DataMat.DisplayUnits = 'mumol.l-1';
    

    all_event_names = {DataMat.Events.label};
    events_found = ismember(selected_event_names, all_event_names);
    if ~all(events_found)
        bst_error(sprintf('Event names "%s" not found (available events: "%s")', ...
                          strjoin(selected_event_names(~events_found), ', '), ...
                          strjoin(all_event_names, ',')));
        return;
    end
    ievents = cellfun(@(l) find(strcmp(l,all_event_names)), selected_event_names);
        
    %% Create model
    
    %Init GLM results struct
    results = []; 
    for data_type = data_types
        bst_progress('text', sprintf('Fitting the model for %s',data_type{1})); 

        % Initialize model 
        model=nst_glm_initialize_model(DataMat.Time);
        
        % Add event-related regressors
        hrf_duration = 32;
        model=nst_glm_add_regressors(model,'event', DataMat.Events(ievents), basis_choice,hrf_duration);
        model=nst_glm_add_regressors(model,'constant');
    
        if  sProcess.options.lfO.Value>=2
            model=nst_glm_add_regressors(model,'linear');
        end
        
        if  sProcess.options.lfO.Value>=3
            lfo_cutoff=sProcess.options.lfo_cutoff.Value{1};
            model=nst_glm_add_regressors(model,'DCT',[1/model.time(end) 1/lfo_cutoff],{'LFO'});  
        end
        
        % Comment this section as DCT is adding too many regressor ( >1500 when
        % using the defaults bands [.2 .6; .8 .15]
        %if sProcess.options.use_DCT.Value 
        %    bands_name=sProcess.options.freqbands.Value(:,1)';
        %    freq_bands=process_tf_bands('GetBounds',sProcess.options.freqbands.Value);
        %    model=nst_glm_add_regressors(model,"DCT",freq_bands,bands_name);  
        %end
        % Add Regressor from file2
        if ~isempty(sInput_ext) && ~isempty(sInput_ext.FileName)
            model=nst_glm_add_regressors(model,'external_input',data_type);
        end
        
        % apply the same filter that was applied to the data to the design
        % matrix
        filter_type = sProcess.options.filter_model.Value;
        low_cutoff = sProcess.options.hpf_low_cutoff.Value{1};
        high_cutoff = sProcess.options.hpf_high_cutoff.Value{1};
        apply_filer = 1; 
        if strcmp(filter_type, 'FIR_bp')
            param = sProcess.options.hpf_transition_band.Value{1};
        elseif strcmp(filter_type, 'IIR_bp')
            param = sProcess.options.hpf_order.Value{1};
        else
            apply_filer = 0;
        end
        if apply_filer
            model = nst_glm_apply_filter(model,filter_type, low_cutoff,high_cutoff,param  );
        end
        
        if sProcess.options.has_detrending.Value
                model = nst_glm_apply_filter(model,'DCT_filter', sProcess.options.dct_cutoff.Value{1}  );
        end

        % Include short-seperation channel
        % Note: we add them after applying filter to the design matrix as
        % we know they have already been filter so we don't filter twice

        if sProcess.options.SS_chan.Value==2 % based on distance
            separation_threshold_m = sProcess.options.SS_chan_distance.Value{1} / 100;
            model=nst_glm_add_regressors(model,'channel',sInput,'distance', separation_threshold_m,data_type);
    
        elseif sProcess.options.SS_chan.Value==3 % based on name  
            if ~isempty(sProcess.options.SS_chan_name.Value)
                SS_name=split(sProcess.options.SS_chan_name.Value,',');
                model=nst_glm_add_regressors(model,'channel',sInput,'name',SS_name',data_type);
            end    
        end   

        % Normalize regressors 
        % model = nst_glm_normalize_design(model);

        % Display Model
        nst_glm_display_model(model,'timecourse');
        
        nb_regressors = size(model.X, 2);

        if isempty(results)
            results = struct('B',zeros(nb_regressors,n_voxel), ...
                             'covB', zeros(nb_regressors,nb_regressors,n_voxel), ...
                             'dfe' , zeros(1,n_voxel) , ...
                             'residuals', zeros( size(model.X,1),n_voxel), ...
                             'mse_residuals', zeros(1,n_voxel) );
        end
        if ~surface_data
            mask    = find(strcmp({ChannelMat.Channel.Group},data_type));
            Y_trim  = Y(:, mask);
        else
            Y_trim = Y;
        end
    
        %% Solve Y = XB + e 
        if sProcess.options.statistical_processing.Value == 1 % Pre-coloring
            method_name = 'OLS_precoloring';
        else
            method_name = 'OLS_prewhitening';
        end
        
        fitted_model= nst_glm_fit(model, Y_trim, method_name);
        results = nst_misc_unpack_glm_result(results, fitted_model,method_name,mask);
    end


    iStudy = sInput.iStudy;
    %% Save results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main output is beta maps stacked in temporal axis
    output_prefix = [sInput.Comment ' | GLM ' method_name ' '];

    sStudy = bst_get('Study', sInput.iStudy);
    
    % Extra fields to save GLM outputs for later internal use
    extra_output.X = model.X; % nb_samples x nb_regressors
    extra_output.X_time = DataMat.Time;
    extra_output.reg_names = model.reg_names;
    extra_output.beta_cov = results.covB; % nb_regressors x nb_regressors x nb_positions
    extra_output.edf = results.dfe;
    extra_output.mse_residuals = results.mse_residuals;
    extra_output.DisplayUnits = DataMat.DisplayUnits; %TODO: check scaling
    
    output_comment = [output_prefix 'fitted model'];
    
    if surface_data
        [sStudy, ResultFile] = nst_bst_add_surf_data(results.B', 1:nb_regressors, [], 'surf_glm_res', output_comment, ...
                                                     [], sStudy, 'GLM', DataMat.SurfaceFile, 0, extra_output);
        OutputFiles{end+1} = ResultFile;
    else
        sDataOut = db_template('data');
        sDataOut.F            = results.B';
        sDataOut.Comment      = output_comment;
        sDataOut.ChannelFlag  = DataMat.ChannelFlag;
        sDataOut.Time         = 1:nb_regressors;
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        
        % Add extra fields
        extra_fields = fieldnames(extra_output);
        for ifield = 1:length(extra_fields)
            sDataOut.(extra_fields{ifield}) = extra_output.(extra_fields{ifield});
        end

        % Save to bst database
        glm_fn = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_glm');
        sDataOut.FileName = file_short(glm_fn);
        bst_save(glm_fn, sDataOut, 'v7');
        % Register in database
        db_add_data(sInput.iStudy, glm_fn, sDataOut);
        OutputFiles{end+1} = glm_fn;
    end
    
    if save_betas
        % Saving B as maps
        for i_reg_name=1:length(model.reg_names)
            data_out = zeros(size(DataMat.F, 1), 1);

            output_tag = sprintf('ir%d_beta%d', sInput.iItem, i_reg_name);
            output_comment = [output_prefix '- beta ' model.reg_names{i_reg_name}];

            if surface_data
                [sStudy, ResultFile] = nst_bst_add_surf_data(results.B(i_reg_name,:)', [1], [], output_tag, output_comment, ...
                                                             [], sStudy, 'GLM', DataMat.SurfaceFile);
            else
                data_out(nirs_ichans,:) = results.B(i_reg_name,:);
                sDataOut = db_template('data');
                sDataOut.F            = data_out;
                sDataOut.Comment      = output_tag;
                sDataOut.ChannelFlag  = DataMat.ChannelFlag;
                sDataOut.Time         = [1];
                sDataOut.DataType     = 'recordings';
                sDataOut.nAvg         = 1;
                sDataOut.DisplayUnits = DataMat.DisplayUnits; %TODO: check scaling

                beta_fn = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['data_beta_' model.reg_names{i_reg_name}]);
                sDataOut.FileName = file_short(beta_fn);
                bst_save(beta_fn, sDataOut, 'v7');
                % Register in database
                db_add_data(sInput.iStudy, beta_fn, sDataOut);
                OutputFiles{end+1} = beta_fn;
            end 
        end
    end
    
    if save_residuals
        % Saving the residual matrix.
        output_tag = sprintf('ir%d_glm_res', sInput.iItem);
        output_comment = [output_prefix '- residuals'];
        if surface_data
                [sStudy, ResultFile] = nst_bst_add_surf_data(results.residuals', DataMat.Time, [], output_tag, output_comment, ...
                                                             [], sStudy, 'GLM', DataMat.SurfaceFile);
                OutputFiles{end+1} = ResultFile;
        else
            
            Out_DataMat = db_template('data');
            Out_DataMat.F           =  results.residuals';
            Out_DataMat.Comment     = output_comment;
            Out_DataMat.DataType     = 'recordings';
            Out_DataMat.Time        =  DataMat.Time;
            Out_DataMat.Events      =  DataMat.Events;
            Out_DataMat.ChannelFlag =  DataMat.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
            Out_DataMat.DisplayUnits = DataMat.DisplayUnits;
            Out_DataMat.nAvg         = 1;
            
            Out_DataMat = bst_history('add', Out_DataMat, DataMat.History, '');
            Out_DataMat = bst_history('add', Out_DataMat, 'GLM', FormatComment(sProcess));
            residuals_fn = bst_process('GetNewFilename', fileparts(sInput.FileName), 'data_residual');
            Out_DataMat.FileName = file_short(residuals_fn);
            bst_save(residuals_fn, Out_DataMat, 'v7');
            db_add_data(iStudy, residuals_fn, Out_DataMat);
            OutputFiles{end+1} = residuals_fn;
            clear residuals;
        end
    end
    
    if save_fit
        fit = model.X*results.B;
        output_tag = sprintf('ir%d_glm_fit', sInput.iItem);
        output_comment = [output_prefix '- signal fit'];
        if surface_data
                [sStudy, ResultFile] = nst_bst_add_surf_data(fit', DataMat.Time, [], output_tag, output_comment, ...
                                                             [], sStudy, 'GLM', DataMat.SurfaceFile);
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
            Out_DataMat = bst_history('add', Out_DataMat, 'GLM', FormatComment(sProcess));
            fit_fn = bst_process('GetNewFilename', fileparts(sInput.FileName), 'data_fit');
            Out_DataMat.FileName = file_short(fit_fn);
            bst_save(fit_fn, Out_DataMat, 'v7');
            db_add_data(iStudy, fit_fn, Out_DataMat);
            OutputFiles{end+1} = fit_fn;
        end        
    end
end


function freqBands= getDefaultFreqBands() 
    freqBands=[ {'heart'},{'.2, .6'},{'all'}; ...
            {'respiration'},{'.8, 1.5'},{'all'} ];
    
end


function lpf = lpf_hrf(h, signal_length)
    % From NIRS_SPM / NIRS10
    h = [h; zeros(size(h))];
    g = abs(fft(h));
    h = real(ifft(g));
    h = fftshift(h)';
    n = length(h);
    d = (1:n) - n/2 -1;
    lpf = spdiags(ones(signal_length,1)*h, d, signal_length, signal_length);
    lpf = spdiags(1./sum(lpf')', 0, signal_length, signal_length) * lpf;
end

function hrf_types = get_hrf_types()
    hrf_types.CANONICAL = 1;
    hrf_types.GAMMA  = 2;
    hrf_types.DECONV = 3;

end

function [X, names, hrf] = make_design_matrix(time, events, hrf_type, hrf_duration, include_trend)
	 warning('Deprecated process');
    if nargin < 4
        hrf_duration = 25;
    end
    
    if nargin < 5
        include_trend = 1;
    end

    X = [];
    names = {};
    
    dt = time(2)-time(1);
    
    %% Setup HRF
    hrf_time = 0:dt:hrf_duration;
    if hrf_time(end) ~= hrf_duration
        warning('HRF duration mismatch due to sampling: %f sec', ...
                hrf_duration-hrf_time(end));
    end
    
    hrf_types = get_hrf_types();
    switch hrf_type
        case hrf_types.CANONICAL  
            hrf = cpt_hrf_canonical(hrf_time); 
        case hrf_types.GAMMA 
            hrf = cpt_hrf_gamma(hrf_time);
        case hrf_types.BOXCAR 
            hrf = cpt_hrf_boxcar(hrf_time);
        otherwise
            bst_error('Unknown hrf_type');
            return;
    end
    
    if size(hrf, 2) ~= 1
        hrf = hrf'; % ensure column vector
    end
    
    %% Make stimulus-induced design matrix
    X = nst_make_event_regressors(events, hrf, time);
    names = {events.label};
    
    %% Add trend function
    if include_trend
        [C,name] = getTrend(time, 'Constant');
        X = [X C];
        names = [names name];
    end
    
    %% Sanity checks
    % Check the rank of the matrix
    if  rank(X) < size(X,2)
        bst_report('Warning', sProcess, sInput, 'The design matrix is not full-ranked');
    end
    
    % Check the collinearity of the matrix
    if cond(X) > 300 
        bst_report('Warning', sProcess, sInput, [ 'The design matrix is high-correlated : Cond(x)=' num2str(cond(X))] );
    end  
    
end

function [C,names]=getTrend(time, trend_choice)
    warning('Deprecated process');
    switch trend_choice
        case 'Constant'
            trend_function=@Constant;
            names={'Constant'};
        otherwise
            trend_function=@Constant;
            names={'Constant'};
    end  
    
    C = trend_function(time);
      
end



function signal = cpt_hrf_gamma(t,peakTime,peakDisp) 
% Gamma : apply the gamma function over t_vect
    
    %signal=zeros( size(t_vect) );  
    if nargin <  3, peakDisp = 1; end
    if nargin <  2, peakTime = 4; end

    signal = peakDisp^peakTime*t.^(peakTime-1).*exp(-peakDisp*t)/gamma(peakTime);
    signal = signal / sum(signal);
end

function signal = cpt_hrf_boxcar(t_vect,lag,duration)
% BoxCar : apply the box car function over t_vect
    
    if nargin <  3, duration = 5; end
    if nargin <  2, lag = 3; end


    signal=zeros( size(t_vect) );
    i=1;
    
    for t=t_vect
        if( t >= lag && t <= lag + duration ) 
            signal(i)=1;
        else
            signal(i)=0;
        end
        i=i+1;
    end
end

function signal= cpt_hrf_canonical_spm(t) 
signal = spm_hrf(diff(t(1:2)));
end

function signal= cpt_hrf_canonical(t,peakTime,uShootTime,peakDisp,uShootDisp,ratio) 
% Canonical :return the Canonical Hrf
    
    assert( isvector(t)  )

    %signal=zeros( size(t_vect) );  
    if nargin <  2, peakTime    = 6.2;end % maybe 6.2 to take into account small delay compared to spm_hrf. TODO: clarify
    if nargin <  3, uShootTime  = 16;end
    if nargin <  4, peakDisp    = 1;end
    if nargin <  5, uShootDisp  = 1;end
    if nargin <  6, ratio       = 1/6;end
          
   signal = peakDisp^peakTime*t.^(peakTime-1).*exp(-peakDisp*t)/gamma(peakTime) - ratio*uShootDisp^uShootTime*t.^(uShootTime-1).*exp(-uShootDisp*t)/gamma(uShootTime);
   signal = signal / sum(signal);
end

function signal= Constant(t,value) 
% Constant : return a constant function.
     warning('Deprecated process'); 
    if nargin <  2, value = 1; end

    assert( isvector(t)  )
    signal = value * ones(numel(t),1);
    
end




