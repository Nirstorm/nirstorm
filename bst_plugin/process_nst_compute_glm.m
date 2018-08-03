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
    sProcess.Comment     = 'Compute GLM';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1401;
    sProcess.isSeparator = 0;

    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-implementation';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','raw'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    sProcess.options.stim_events.Comment = 'Stim events: ';
    sProcess.options.stim_events.Type    = 'text';
    sProcess.options.stim_events.Value   = '';
    
    sProcess.options.basis.Comment = 'Basis function';
    sProcess.options.basis.Type    = 'combobox';
    sProcess.options.basis.Value   = {1, {'Hrf','Gamma','BoxCar'}};
    
    sProcess.options.trend.Comment = 'Add a constant trend function ';
    sProcess.options.trend.Type    = 'checkbox';
    sProcess.options.trend.Value   =  1;
    
    sProcess.options.external.Comment = 'Call external function to define external regressor';
    sProcess.options.external.Type    = 'text';
    sProcess.options.external.Value   =  '';    
    
    
    sProcess.options.save.Comment = 'Save the design matrix';
    sProcess.options.save.Type    = 'checkbox';
    sProcess.options.save.Value   =  1;
    
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    sProcess.options.fitting.Comment = 'Fitting method';
    sProcess.options.fitting.Type    = 'combobox';
    sProcess.options.fitting.Value   = {1, {'OLS','AR-IRLS(AnalyzIR toolbox)' }};

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
    fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
    if( fitting_choice == 1 )
        Comment=[ Comment ' OLS fit'];
    elseif( fitting_choice == 2)
        Comment=[ Comment ' AR-IRLS fit'];
    end    
end

function OutputFiles = Run(sProcess, sInput)

    OutputFiles={};

    DataMat = in_bst_data(sInput.FileName);
    basis_choice=sProcess.options.basis.Value{2}(sProcess.options.basis.Value{1});

    %% Select events
    % TODO: utests with all events found, no event selected, no available
    %       event in data, one event missing
    if isempty(sProcess.options.stim_events.Value)
         bst_error('No event selected');
    end
    selected_event_names = cellfun(@strtrim, strsplit(sProcess.options.stim_events.Value, ','),...
                                   'UniformOutput', 0);
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
    [X,names]=getX(DataMat.Time, DataMat.Events(ievents), basis_choice);
    
    % Add trend function
    if sProcess.options.trend.Value == 1 
        [C,name]=getTrend(DataMat.Time, 'Constant');
        X=[X C];
        names=[names name];
    end
    if isfield(sProcess.options,'external') && ~strcmp(sProcess.options.external.Value,'')
        if strcmp(sProcess.options.external.Value(1),'@') % this is a function 
            external_function_names = sProcess.options.external.Value;
           try
               [C,name]=feval( str2func(external_function_names),sProcess, sInput);

 
           catch
               bst_error([ 'Error during the call of ' external_function_names ]);  
               return;
           end
        else % this is a .mat  file
            if exist(sProcess.options.external.Value,'file')
                data=load(sProcess.options.external.Value);
                if isfield(data,'C') &&  isfield(data,'name')
                    C=data.C;
                    name=data.name;
                else
                    bst_error([ 'Can''t find structues C and names in ' external_function_names ]);  
                    return;
                end
                
            else
                bst_error([ external_function_names 'not found']);  
                return;    
            end
            
        end
        
        if ~(size(C,1) == size(X,1))
            bst_error([ 'Dimension of the external regressor returned by ' external_function_names ...
                         'doesn''t match design matrix dimension']);
            return;
        end
        if ~(size(C,2) == length(name))
            bst_error([external_function_names ' have to return one name for each regressor']);    
            return;
        end
        X=[X C];
        names=[names name];             
    end
    
    % Check the rank of the matrix
    n_regressor=size(X,2);
    if  rank(X) < n_regressor
        bst_report('Warning', sProcess, sInput, 'The design matrix is not full-ranked');
    end
    
    % Check the collinearity of the matrix
    
    if cond(X) > 300 
        bst_report('Warning', sProcess, sInput, [ 'The design matrix is high-correlated : Cond(x)=' num2str(cond(X))] );
    end    

    
    iStudy = sInput.iStudy;

    if sProcess.options.save.Value == 1 
        % Save the design matrix
        Out_DataMat = db_template('matrixmat');
        Out_DataMat.Value           = X';
        Out_DataMat.Comment     = 'Design Matrix';
        Out_DataMat.Description = names'; % Names of the regressors 
        Out_DataMat.ChannelFlag =  ones(n_regressor,1);   % List of good/bad channels (1=good, -1=bad)
        Out_DataMat.Time        =  DataMat.Time;
        Out_DataMat.DisplayUnits = DataMat.DisplayUnits; 
        Out_DataMat = bst_history('add', Out_DataMat, 'Design', FormatComment(sProcess));
        
        OutputFiles{end+1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'design_matrix');
        bst_save(OutputFiles{end}, Out_DataMat, 'v6');
        db_add_data(iStudy, OutputFiles{end}, Out_DataMat);
    end 
    
    %% Solving  B such as Y = XB +e 
    
    % Get signals of NIRS channels only:
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');
    Y = DataMat.F(nirs_ichans,:)';

    fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
    if( fitting_choice == 1 ) % Use OLS : : \( B= ( X^{T}X)^{-1} X^{T} Y \)
         method_name='OLS_fit';
         [B,covB,dfe]=nst_ols_fit( Y, X );
    elseif( fitting_choice == 2 ) % Use AR-IRLS 
        analyzIR_url='https://bitbucket.org/huppertt/nirs-toolbox/';
        
        if ~exist('nirs.core.ChannelStats','class')
        	bst_error(['AnalyzIR toolbox required. See ' analyzIR_url]);
            return
        end
        method_name='AR-IRLS_fit';
        [B,covB,dfe]=nst_ar_irls_fit( Y, X, round(4/(DataMat.Time(2)-DataMat.Time(1))) ); 
    else
        bst_report('Error', sProcess, sInputs, 'This method is not implemented yet' );
    end
    
    residual=Y - X*B ;
    
   
    %% Save results
    
    % Saving the B Matrix
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value           = B';
    Out_DataMat.Comment     = ['GLM_beta_matrix_' method_name];
    Out_DataMat.Description = names; % Names of the regressors 
    Out_DataMat.ChannelFlag =  ones(size(B,2),1);   % List of good/bad channels (1=good, -1=bad)
    Out_DataMat = bst_history('add', Out_DataMat, DataMat.History, '');
    Out_DataMat = bst_history('add', Out_DataMat, 'GLM', FormatComment(sProcess));

    OutputFiles{end+1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'B_matrix');
    bst_save(OutputFiles{end}, Out_DataMat, 'v6');
    db_add_data(iStudy, OutputFiles{end}, Out_DataMat);
    
    
    % Saving B as maps
    for i_reg_name=1:length(names)
        data_out = zeros(size(DataMat.F, 1), 1);
        data_out(nirs_ichans,:) = B(i_reg_name,:);
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = ['GLM_beta_' method_name '_' names{i_reg_name}];
        sDataOut.ChannelFlag  = DataMat.ChannelFlag;
        sDataOut.Time         = [1];
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = 'U.A.'; %TODO: check scaling
        
        OutputFiles{end+1} = bst_process('GetNewFilename', fileparts(sInput.FileName), ['data_beta_' names{i_reg_name}]);
        bst_save(OutputFiles{end}, sDataOut, 'v6');
        db_add_data(iStudy, OutputFiles{end}, sDataOut);
    end
    
    % Saving the covB Matrix.
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value           = covB;
    Out_DataMat.Comment     = [ 'GLM_covB_' method_name '_df=' int2str(dfe) ];
    Out_DataMat = bst_history('add', Out_DataMat, DataMat.History, '');
    Out_DataMat = bst_history('add', Out_DataMat, 'GLM', FormatComment(sProcess));
    
    OutputFiles{end+1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'covB_matrix');
    bst_save(OutputFiles{end}, Out_DataMat, 'v6');
    db_add_data(iStudy, OutputFiles{end}, Out_DataMat);
        
   
    
    % Saving the residual matrix.
    Out_DataMat = db_template('data');
    Out_DataMat.F           = residual' ;
    Out_DataMat.Comment     = ['GLM_' method_name '_Residual Matrix'];
    Out_DataMat.DataType     = 'recordings'; 
    Out_DataMat.Time        =  DataMat.Time;
    Out_DataMat.Events      =  DataMat.Events;
    Out_DataMat.ChannelFlag =  DataMat.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
    Out_DataMat.DisplayUnits = DataMat.DisplayUnits; 
    Out_DataMat.nAvg         = 1; 
    
    Out_DataMat = bst_history('add', Out_DataMat, DataMat.History, '');
    Out_DataMat = bst_history('add', Out_DataMat, 'GLM', FormatComment(sProcess));
    
    OutputFiles{end+1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'data_residual');
    bst_save(OutputFiles{end}, Out_DataMat, 'v6');
    db_add_data(iStudy, OutputFiles{end}, Out_DataMat);
    

end


function [X,names]=getX(time,events,basis_choice)
	n_event=length(events);
    n_sample=length(time);
    
    X=zeros(n_sample,n_event); 
    
    % removing offset
    time_offset=time(1);   
    sample_offset=round(time(1)/(time(2)-time(1)));
    
    time=time-time_offset;

    % Selecting the basis function 
    switch cell2mat(basis_choice)
        case 'Hrf'  
            basis_function=@Canonical; 
        case'Gamma' 
            basis_function=@Gamma;
        case 'BoxCar' 
            basis_function=@BoxCar;
        otherwise
            basis_function=@Canonical;
    end
    
    % create the stim matrix    
    for i=1:n_event
       %disp(events(i).label)
       names{i}=events(i).label;
       n_run=length(events(i).samples);
       for j =  1:n_run
           event=(events(i).samples(1,j) - sample_offset ):(events(i).samples(2,j) - sample_offset);
           X(event,i)=ones(size(event)); 
       end
    end
    
    
    for i=1:n_event
        X(:,i) = filter(basis_function(time ), 1, X(:,i));
        %X(:,i) = X(:,i) - mean(X(:,i));
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



function signal= Gamma(t,peakTime,peakDisp) 
% Gamma : apply the gamma function over t_vect
    
    %signal=zeros( size(t_vect) );  
    if nargin <  3, peakDisp = 1; end
    if nargin <  2, peakTime = 4; end

    signal = peakDisp^peakTime*t.^(peakTime-1).*exp(-peakDisp*t)/gamma(peakTime);
    signal = signal / sum(signal);
end

function signal= BoxCar(t_vect,lag,duration)
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

function signal= Canonical(t,peakTime,uShootTime,peakDisp,uShootDisp,ratio) 
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




