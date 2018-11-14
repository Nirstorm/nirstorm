function varargout = process_nst_compute_glm( varargin )
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
% Authors: Edouard Delaire, Thomas Vincent 2018
%  
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - design and fit';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1401;
    sProcess.isSeparator = 0;

    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw', 'results'};
    sProcess.OutputTypes = {'data', 'data', 'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.stim_events.Comment = 'Stimulation events: ';
    sProcess.options.stim_events.Type    = 'text';
    sProcess.options.stim_events.Value   = '';
    
    sProcess.options.hrf_model.Comment = 'HRF model: ';
    sProcess.options.hrf_model.Type    = 'combobox';
    sProcess.options.hrf_model.Value   = {1, fieldnames(get_hrf_types())};
    
    sProcess.options.trend.Comment = 'Add constant regressor ';
    sProcess.options.trend.Type    = 'checkbox';
    sProcess.options.trend.Value   =  1;
    
    sProcess.options.fitting.Comment = 'Fitting method';
    sProcess.options.fitting.Type    = 'combobox';
    sProcess.options.fitting.Value   = {1, {'OLS','AR-IRLS' }};

    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    sProcess.options.save_residuals.Comment = '<B>Extra outputs</B>:';
    sProcess.options.save_residuals.Type    = 'label';
    
    sProcess.options.save_residuals.Comment = 'Residuals';
    sProcess.options.save_residuals.Type    = 'checkbox';
    sProcess.options.save_residuals.Value   =  0;
    
    sProcess.options.save_betas.Comment = 'Beta maps';
    sProcess.options.save_betas.Type    = 'checkbox';
    sProcess.options.save_betas.Value   =  0;

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
%     fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
%     if( fitting_choice == 1 )
%         Comment=[ Comment ' OLS fit'];
%     elseif( fitting_choice == 2)
%         Comment=[ Comment ' AR-IRLS fit'];
%     end    
end

function OutputFiles = Run(sProcess, sInput)

    OutputFiles={};

    DataMat = in_bst_data(sInput.FileName);
    basis_choice = sProcess.options.hrf_model.Value{1};

    save_residuals = sProcess.options.save_residuals.Value;
    save_betas = sProcess.options.save_residuals.Value;
    
    %% Select events
    % TODO: utests with all events found, no event selected, no available
    %       event in data, one event missing
    if isempty(sProcess.options.stim_events.Value)
         bst_error('No event selected');
    end
    selected_event_names = cellfun(@strtrim, strsplit(sProcess.options.stim_events.Value, ','),...
                                   'UniformOutput', 0);
                               
    if isfield(DataMat, 'SurfaceFile')
        surface_data = 1;
        parent_data = in_bst_data(DataMat.DataFile);
        % Make sure time axis is consistent
        assert(all(parent_data.Time == DataMat.Time));
        if isempty(DataMat.Events) && isfield(parent_data, 'Events')
            DataMat.Events = parent_data.Events;
        end
        if isempty(DataMat.F) && ~isempty(DataMat.ImageGridAmp) && size(DataMat.ImageGridAmp, 2)==length(DataMat.Time)
            Y = DataMat.ImageGridAmp'; %TODO: check that it doesn't take too much memory
        else
            bst_error('Cannot get signals from surface data');
        end
    else
        surface_data = 0;

        % Get signals of NIRS channels only:
        ChannelMat = in_bst_channel(sInput.ChannelFile);
        [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');
        Y = DataMat.F(nirs_ichans,:)';
    end
    
    all_event_names = {DataMat.Events.label};
    events_found = ismember(selected_event_names, all_event_names);
    if ~all(events_found)
        bst_error(sprintf('Event names "%s" not found in "%s"', ...
                          strjoin(selected_event_names(~events_found), ', '), ...
                          strjoin(all_event_names, ',')));
        return;
    end
    ievents = cellfun(@(l) find(strcmp(l,all_event_names)), selected_event_names);
    
    % Check that all selected events are extended
    % TODO: also handle single events for event-related design
    isExtended = false(1,length(ievents));
    for ievt=1:length(ievents)
        isExtended(ievt) = (size(DataMat.Events(ievents(ievt)).samples, 1) == 2);
    end
    if ~all(isExtended)
         bst_error(sprintf('Simple events not supported: %s ', ...
                           strjoin(selected_event_names(ievents(~isExtended)), ', ')));
         return;
    end
    
    %% Create the design matrix X
    hrf_duration = 25; % sec -- TODO: expose as process parameter
    [X,names] = make_design_matrix(DataMat.Time, DataMat.Events(ievents), ...
                                   basis_choice, hrf_duration, sProcess.options.trend.Value);  
    n_regressor = size(X, 2);
    
    iStudy = sInput.iStudy;

    
    
    %% Solving  B such as Y = XB +e 

    fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
    if( fitting_choice == 1 ) % Use OLS : : \( B= ( X^{T}X)^{-1} X^{T} Y \)
         method_name='OLS_fit';
         [B,covB,dfe]=ols_fit( Y, X );
    elseif( fitting_choice == 2 ) % Use AR-IRLS 
        analyzIR_url='https://bitbucket.org/huppertt/nirs-toolbox/';
        
        if ~exist('nirs.core.ChannelStats','class')
        	bst_error(['AnalyzIR toolbox required. See ' analyzIR_url]);
            return
        end
        method_name='AR-IRLS_fit';
        [B,covB,dfe]=ar_irls_fit( Y, X, round(4/(DataMat.Time(2)-DataMat.Time(1))) ); 
    else
        bst_report('Error', sProcess, sInput, 'This method is not implemented yet' );
    end
    
    output_prefix = [sInput.Comment ' | GLM ' method_name ' '];
    
    residual = Y - X*B ;
    
   
    %% Save results
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value       = X'; % nb_regressors x nb_samples
    Out_DataMat.Comment     = [output_prefix 'results'];
    Out_DataMat.Description = names'; % Names of the regressors
    Out_DataMat.ChannelFlag = DataMat.ChannelFlag;  
    Out_DataMat.Time        = DataMat.Time;
    Out_DataMat.DataFile    = sInput.FileName;
    Out_DataMat.SurfaceFile = DataMat.SurfaceFile;
    
    Out_DataMat.beta = B; % nb_regressors x nb_positions
    Out_DataMat.beta_cov = covB; % nb_regressors x nb_regressors x nb_positions
    Out_DataMat.df = dfe; % scalar int
    Out_DataMat.input_units = DataMat.DisplayUnits;
    
    Out_DataMat = bst_history('add', Out_DataMat, 'Design', FormatComment(sProcess));
    OutputFiles{end+1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'design_matrix');
    % Save on disk
    save(OutputFiles{end}, '-struct', 'Out_DataMat');
    % Register in database
    db_add_data(iStudy, OutputFiles{end}, Out_DataMat);
    
    if save_betas
        % Saving B as maps
        for i_reg_name=1:length(names)
            data_out = zeros(size(DataMat.F, 1), 1);

            output_tag = sprintf('ir%d_beta%d', sInput.iItem, i_reg_name);
            output_comment = [output_prefix '- beta ' names{i_reg_name}];
            sStudy = bst_get('Study', sInput.iStudy);

            if surface_data
                [sStudy, ResultFile] = nst_bst_add_surf_data(B(i_reg_name,:)', [1], [], output_tag, output_comment, ...
                                                             sInput, sStudy, 'GLM', DataMat.SurfaceFile);
            else
                data_out(nirs_ichans,:) = B(i_reg_name,:);
                sDataOut = db_template('data');
                sDataOut.F            = data_out;
                sDataOut.Comment      = output_tag;
                sDataOut.ChannelFlag  = DataMat.ChannelFlag;
                sDataOut.Time         = [1];
                sDataOut.DataType     = 'recordings';
                sDataOut.nAvg         = 1;
                sDataOut.DisplayUnits = sInput.DisplayUnits; %TODO: check scaling

                beta_fn = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['data_beta_' names{i_reg_name}]);
                sDataOut.FileName = file_short(beta_fn);
                bst_save(beta_fn, sDataOut, 'v7');
                % Register in database
                db_add_data(sInput.iStudy, beta_fn, sDataOut);
            end 
        end
    end
    
    if save_residuals
        % Saving the residual matrix.
        output_tag = sprintf('ir%d_glm_res', sInput.iItem);
        output_comment = [output_prefix '- residuals'];
        if surface_data
                [sStudy, ResultFile] = nst_bst_add_surf_data(residual', DataMat.Time, [], output_tag, output_comment, ...
                                                             sInput, sStudy, 'GLM', DataMat.SurfaceFile);
        else
            Out_DataMat = db_template('data');
            Out_DataMat.F           = residual' ;
            Out_DataMat.Comment     = output_tag;
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
        end
    end
end


function [B,covB,dfe]=ols_fit(y,X)
    
    B= pinv(X'*X)* X'*y; % or B=X\y; 
    
    dfe = size(y,1) - rank(X);
    
    fit = X * B;
    residual = y - fit;     
    
    % covE=cov(residual);
    
    mse_residuals = var(residual) * (size(y,1)-1) / dfe;
    
    inv_XtX = pinv(transpose(X)*X);
    for i=1:size(residual,2)
        covB(:,:,i) = mse_residuals(i) * inv_XtX;
        % covB(:,:,i) = covE(i,i) * pinv(transpose(X)*X);
    end
    

    
    if 0 % for test when running GLMTest.test_cortical_simulation
       
        % activ pos hbO: 4708
        % inactiv pos hbo: 4977
        
        poi_inact = 4955;
        
        figure(); hold on;
        plot(y(:,poi_inact), 'k');
        plot(residual(:,poi_inact), 'g');
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

function [B,covB,dfe]=ar_irls_fit(y,X,pmax)
	stat=nirs.math.ar_irls(y,X, pmax );
    
    B=stat.beta;
    
    covB=zeros( size(stat.covb,1), size(stat.covb,2),size(stat.covb,3));
    for i=1:size(stat.covb,3)
       covB(:,:,i)= stat.covb(:,:,i,i);
    end    
    dfe=stat.dfe;
end

function hrf_types = get_hrf_types()
hrf_types.CANONICAL = 1;
hrf_types.GAMMA = 2;
hrf_types.BOXCAR = 3;
end

function [X,names] = make_design_matrix(time, events, hrf_type, hrf_duration, include_trend)
	
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
    n_samples = length(time);
    X = nst_make_event_regressors(events, hrf, n_samples);
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
        end;
        i=i+1;
    end
end

function signal= cpt_hrf_canonical(t,peakTime,uShootTime,peakDisp,uShootDisp,ratio) 
% Canonical :return the Canonical Hrf
    
    assert( isvector(t)  )

    %signal=zeros( size(t_vect) );  
    if nargin <  2, peakTime    = 4;end
    if nargin <  3, uShootTime  = 16;end
    if nargin <  4, peakDisp    = 1;end
    if nargin <  5, uShootDisp  = 1;end
    if nargin <  6, ratio       = 1/6;end
          
   signal = peakDisp^peakTime*t.^(peakTime-1).*exp(-peakDisp*t)/gamma(peakTime) - ratio*uShootDisp^uShootTime*t.^(uShootTime-1).*exp(-uShootDisp*t)/gamma(uShootTime);
   signal = signal / sum(signal);
end

function signal= Constant(t,value) 
% Constant : return a constant function.
    
    %signal=zeros( size(t_vect) );  
    if nargin <  2, value = 1; end

    assert( isvector(t)  )
    signal = value * ones(numel(t),1);
    
end




