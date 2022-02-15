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
    
    sProcess.options.filter_model.Type    = 'radio_line';
    sProcess.options.filter_model.Comment   = {'No filter','IIR filter','Detrending','Process used to remove slow fluctuations: '};
    sProcess.options.filter_model.Value=1;
    
    sProcess.options.hpf_low_cutoff.Comment = 'High-pass filter frequency: ';
    sProcess.options.hpf_low_cutoff.Type    = 'value';
    sProcess.options.hpf_low_cutoff.Value   = {0.01, 'Hz', 2};
    
    sProcess.options.dct_cutoff.Comment = 'Detrend miminum Period (0= only linear detrend): ';
    sProcess.options.dct_cutoff.Type    = 'value';
    sProcess.options.dct_cutoff.Value   = {200, 's', 0};
    
    sProcess.options.trim_start.Comment = 'Ignore starting signal: ';
    sProcess.options.trim_start.Type    = 'value';
    sProcess.options.trim_start.Value   = {0, 'sec', 2};
    sProcess.options.trim_start.Hidden   = 1; % Hide, as it is deprecated
    
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

    DataMat = in_bst_data(sInput.FileName);
    basis_choice = sProcess.options.hrf_model.Value{1};

    save_residuals = sProcess.options.extra_output.Value || sProcess.options.save_residuals.Value;
    save_betas = sProcess.options.extra_output.Value ||  sProcess.options.save_betas.Value;
    save_fit = sProcess.options.extra_output.Value || sProcess.options.save_fit.Value;

    trim_start = sProcess.options.trim_start.Value{1};
    if trim_start > 0
       warning(['trim_start option is deprecated. ' ...
                'Use process_extract_time beforehand instead.']); 
    end
    
    %% Select events
    % TODO: utests with all events found, no event selected, no available
    %       event in data, one event missing
    if isempty(sProcess.options.stim_events.Value)
         bst_error('No event selected');
    end
    selected_event_names = cellfun(@strtrim, strsplit(sProcess.options.stim_events.Value, ','),...
                                   'UniformOutput', 0);
    
    %% Load data and events
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    if isfield(DataMat, 'SurfaceFile')
        surface_data = 1;
        channel_data = in_bst_data(DataMat.DataFile);
        % Make sure time axis is consistent
        assert(all(channel_data.Time == DataMat.Time));
        if isempty(DataMat.Events) && isfield(channel_data, 'Events')
            DataMat.Events = channel_data.Events;
        end
        if isempty(DataMat.F) && ~isempty(DataMat.ImageGridAmp) && size(DataMat.ImageGridAmp, 2)==length(DataMat.Time)
            Y = DataMat.ImageGridAmp';
            if issparse(Y)
                Y = full(Y);
            end
            n_voxel=size(Y,2);
            mask=find(~all(Y==0,1)); % Only keep channel in the field of view
            Y=Y(:,mask);
        else
            bst_error('Cannot get signals from surface data');
        end
    else
        surface_data = 0;
        % Get signals of NIRS channels only:
        [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');
        Y = DataMat.F(nirs_ichans,:)';
        mask=[];
        n_voxel=size(Y,2);
    end
    
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
    
    %% previus processing 
    filter_type='none';
    filter_param=0;
    if sProcess.options.statistical_processing.Value == 2
        filter_type='IIR_highpass';
        filter_param=sProcess.options.hpf_low_cutoff.Value{1};
    elseif sProcess.options.statistical_processing.Value == 3
        filter_type='DCT_filter';
        filter_param=sProcess.options.dct_cutoff.Value{1};  
    end        
    
    
    %% Create model
    hrf_duration = 32; % sec -- TODO: expose as process parameter?
    % Initialize model 
    model=nst_glm_initialize_model(DataMat.Time);
    
    % Add event-related regressors
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
    
    
    % Include short-seperation channel
    if sProcess.options.SS_chan.Value==2 % based on distance
        separation_threshold_m = sProcess.options.SS_chan_distance.Value{1} / 100;
        model=nst_glm_add_regressors(model,'channel',sInput,'distance', separation_threshold_m,{'HbO'});
        model=nst_glm_add_regressors(model,'channel',sInput,'distance', separation_threshold_m,{'HbR'});

    elseif sProcess.options.SS_chan.Value==3 % based on name  
        if ~isempty(sProcess.options.SS_chan_name.Value)
            SS_name=split(sProcess.options.SS_chan_name.Value,',');
            model=nst_glm_add_regressors(model,'channel',sInput,'name',SS_name',{'HbO'});
            model=nst_glm_add_regressors(model,'channel',sInput,'name',SS_name',{'HbR'});
        end    
    end   
    
    % Add Regressor from file2
    if ~isempty(sInput_ext) && ~isempty(sInput_ext.FileName)
                
        hb_types = {'HbO', 'HbR', 'HbT'};
        hb_type = [];
        for ihb=1:length(hb_types)
            if ~isempty(strfind(lower(sInput.FileName), lower(hb_types{ihb})))
                hb_type = hb_types{ihb};
                break;
            end
        end
        model=nst_glm_add_regressors(model,'external_input',hb_type);
    end
    
    % Display Model
    nst_glm_display_model(model,'timecourse');
    
    nb_regressors = size(model.X, 2);
    iStudy = sInput.iStudy;

    dt = diff(DataMat.Time(1:2));
    trim_start_sample = round(trim_start / dt);
    Y_trim = Y((trim_start_sample+1):end, :);
    model.X = model.X((trim_start_sample+1):end, :);

    %% Solve Y = XB + e 
    if sProcess.options.statistical_processing.Value == 1 % Pre-coloring
        method_name = 'OLS_precoloring';
    else
        method_name = 'OLS_prewhitening';
    end
    
    fitted_model= nst_glm_fit(model, Y_trim, filter_type, filter_param, method_name);
    [B,covB,dfe,residuals,mse_residuals] = nst_misc_unpack_glm_result(fitted_model,method_name,surface_data,nb_regressors,n_voxel,mask);
    

    
    %% Save results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main output is beta maps stacked in temporal axis
    output_prefix = [sInput.Comment ' | GLM ' method_name ' '];

    sStudy = bst_get('Study', sInput.iStudy);
    
    % Extra fields to save GLM outputs for later internal use
    extra_output.X = model.X; % nb_samples x nb_regressors
    extra_output.X_time = DataMat.Time;
    extra_output.reg_names = model.reg_names;
    extra_output.beta_cov = covB; % nb_regressors x nb_regressors x nb_positions
    extra_output.edf = dfe;
    extra_output.mse_residuals = mse_residuals;
    extra_output.DisplayUnits = DataMat.DisplayUnits; %TODO: check scaling
    
    output_comment = [output_prefix 'fitted model'];
    
    if surface_data
        [sStudy, ResultFile] = nst_bst_add_surf_data(B', 1:nb_regressors, [], 'surf_glm_res', output_comment, ...
                                                     [], sStudy, 'GLM', DataMat.SurfaceFile, 0, extra_output);
        OutputFiles{end+1} = ResultFile;
    else
        sDataOut = db_template('data');
        sDataOut.F            = B';
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
                [sStudy, ResultFile] = nst_bst_add_surf_data(B(i_reg_name,:)', [1], [], output_tag, output_comment, ...
                                                             [], sStudy, 'GLM', DataMat.SurfaceFile);
            else
                data_out(nirs_ichans,:) = B(i_reg_name,:);
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
                [sStudy, ResultFile] = nst_bst_add_surf_data(residuals', DataMat.Time((trim_start_sample+1):end), [], output_tag, output_comment, ...
                                                             [], sStudy, 'GLM', DataMat.SurfaceFile);
                OutputFiles{end+1} = ResultFile;
        else
            
            Out_DataMat = db_template('data');
            Out_DataMat.F           =  residuals';
            Out_DataMat.Comment     = output_comment;
            Out_DataMat.DataType     = 'recordings';
            Out_DataMat.Time        =  DataMat.Time((trim_start_sample+1):end);
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
        fit = model.X*B;
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

function [B,proj_X] = compute_B_svd(y,X)
    warning('Deprecated process');
    
    [X_u, X_s, X_v] = svd(X,0);
    X_diag_s = diag(X_s);
    X_rank_tol =  max(size(X)) * max(abs(X_diag_s)) * eps;
    X_rank =  sum(X_diag_s > X_rank_tol);

    % Compute projector based on X to solve Y = X*B
    proj_X = X_v(:,1:X_rank) * diag(1./X_diag_s(1:X_rank)) * X_u(:,1:X_rank)';
    
    % Fit to data
    B = proj_X * y; 

end
function [B, covB, dfe, residuals, mse_residuals] = ols_fit(y, dt, X, hrf, hpf_low_cutoff)
    warning('Deprecated process');
    % Low pass filter, applied to input data
    lpf = full(lpf_hrf(hrf, size(y, 1)));
    % Convert to full mat since operation on sparse matrices are
    % single-threaded only but full 
    y_filtered = lpf * y;
    
    if any(~isfinite(y_filtered))
        warning('Found non-finite values in filtered input data.');
    end

    % Band-pass filtering of the design matrix
    if hpf_low_cutoff > 0 
        X_filtered = lpf * process_nst_iir_filter('Compute', X(:,1:(end-1)), 1/dt, ...
                                              'highpass', hpf_low_cutoff, ...
                                               0, 2, 0);
    else
        X_filtered = lpf * X(:,1:(end-1));
    end    
    X_filtered = [X_filtered X(:,end)];
    
    %% Solve y = X*B using SVD as in SPM
    [B,proj_X] = compute_B_svd(y_filtered,X_filtered);
    
    %% For stat afterwards
    res_form_mat = eye(size(lpf)) - X_filtered * proj_X;
    lpf_lpf_T = full(lpf * lpf');
    RV = res_form_mat * lpf_lpf_T;
    trRV = sum(diag(RV));
    RVRVt = RV .* RV';
    trRVRV = sum(RVRVt(:)); % faster than sum(diag(RV * RV));
    dfe = trRV^2 / trRVRV;
    
    fit = X_filtered * B;
    residuals = y_filtered - fit;     
    mse_residuals = var(residuals) * (size(y,1)-1) / trRV;
        
%     figure();plot(y_filtered(:, 161)/max(y_filtered(:,161)), 'b', 'LineWidth', 2); hold on; plot(X_filtered);
%     figure(); hold on; plot(y_filtered(:,161), 'b'); plot(residual(:,161), 'g'); plot(fit(:,161), 'r');
%   
    pXS = proj_X * lpf;
    covB = pXS * pXS';

    if 0 % for test when running GLMTest.test_cortical_simulation
       
        % activ pos hbO: 4708
        % inactiv pos hbo: 4977
        
        poi_inact = 4955;
        
        figure(); hold on;
        plot(y(:,poi_inact), 'k');
        plot(residuals(:,poi_inact), 'g');
        plot(fit(:,poi_inact), 'r');
        mse_residuals_inact = mse_residuals(poi_inact);
        fprintf('MSE_inact = %e\n', mse_residuals_inact);
        t_stat_inact = B(poi_inact) / sqrt(covB(1,1,poi_inact));
        fprintf('tstat_inact = %1.3f\n', t_stat_inact);
        p_val_inact = process_test_parametric2('ComputePvalues', t_stat_inact, dfe, 't', ...
                                               'one+');
        fprintf('p_val_inact = %1.3f\n', p_val_inact);
        
    end
end
function [B_out, covB_out, dfe_out, residuals_out, mse_residuals_out] = AR1_ols_fit(Y, dt, X, hpf_low_cutoff)
    
    warning('Deprecated process');

    n_chan=size(Y,2);
    n_cond=size(X,2);
    n_time=size(Y,1);
    
    max_iter=10;
    B_out=zeros(n_cond,n_chan);
    
    covB_out=zeros(n_cond,n_cond,n_chan);
    dfe_out=zeros(1,n_chan);
    residuals_out=zeros(n_time,n_chan);
    mse_residuals_out=zeros(1,n_chan);
    
    % high-pass filtering of the design matrix
    if hpf_low_cutoff > 0 
        X_hpf = process_nst_iir_filter('Compute', X(:,1:(end-1)), 1/dt, ...
                                              'highpass', hpf_low_cutoff, ...
                                               0, 2, 0);    
    else
        X_hpf = X(:,1:(end-1));
    end
    
    X_hpf = [X_hpf X(:,end)];
    
    %X_hpf=X;
    
    [B_init,proj_X] = compute_B_svd(Y,X_hpf);
    bst_progress('start', 'GLM - Pre-whitenning ' , 'Fitting the GLM', 1, n_chan);
    
    parfor i_chan=1:n_chan
        % Solve B for chan i_chan. We need to solve B for each channel
        % speratly as we are fitting one AR model per channel. This might 
        % not be a good idead in the source space. 
        tic;
        iter=0;
       
        y=Y(:,i_chan);
        SX=X_hpf;
        SY=y;
       
        B=B_init(:,i_chan);
        B0=zeros(n_cond,1);
        S_prod=speye(n_time);
        
        while( norm(B-B0)/norm(B) > 1e-2 && iter < max_iter )
            
            B0=B;
            
            % Estimate the AR(1) processe on the residual
            res=SY-SX*B0; 
            W=nst_math_fit_AR(res',1);
        
            %W=[1  -W(2:end)];
            %X_filtered(:,1:end-1)=filter(1,W,X_filtered(:,1:end-1));
            %y_filtered=filter(1,W,y_filtered);
       
            % Compute the filtering matrix 
            a1 = -W(2);
            ka = sparse( ( eye(n_time) - diag( ones(1,n_time-1)*a1,-1) ))^(-1); 
            Va = ka* ka';
        
            S= full( inv(Va) )^(0.5); % can't use power .5 on sparse matrix
            S_prod=sparse(S)*S_prod;
            
            % Apply the filtering to the data and the deisgn matrix
            SY = S*SY;
            SX = S*SX;

            % Compute B for Sy = SXB + Se following an iid normal distribution
            [B,proj_X] = compute_B_svd(SY,SX);
       
            iter=iter+1;
            disp([ 'iter ' num2str(iter) ' norm(B-B0)/norm(B) = ' num2str(norm(B-B0)/norm(B)) ])
       end
       
        S=S_prod;
        Va=full(inv(S))^2;
        S=full(S);
        
        % Compute stat afterwards
       
        pSX=pinv(SX);
        
        R = eye(n_time) - SX * pSX;
        
        S_Va_S_t = S * Va * S';
        RV = R * S_Va_S_t;
        trRV = sum(diag(RV));
        RVRVt = RV .* RV';
        trRVRV = sum(RVRVt(:)); % faster than sum(diag(RV * RV));
        
        fit = SX * B;
        residuals = SY - fit;
        sigma2=var(residuals);

        dfe = trRV^2 / trRVRV;
        mse_residuals =  sigma2* (size(y,1)-1)/trRV;
        covB = pSX * pSX';

        
        % Save stats 
        
               
        B_out(:,i_chan)=B;
        covB_out(:,:,i_chan)=covB;
        dfe_out(i_chan)=dfe;
        residuals_out(:,i_chan)=residuals;
        mse_residuals_out(i_chan)=mse_residuals;
        e = toc;
        disp( [ '#' num2str(i_chan) ' analized in ' num2str(iter) ' iteration (' num2str(e) ' sec)']) 
        bst_progress('inc', 1); 
    end
    
    bst_progress('stop');

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
hrf_types.GAUSSIAN = 2;
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




