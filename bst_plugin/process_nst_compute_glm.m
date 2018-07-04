function varargout = process_nst_compute_glm( varargin )
% process_compute_glm: compute the glm : find B such as Y = XB +e with X
% the design matrix of the experimentation. 
% This process uses OLS to fit the data.
%
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
% Authors: Edouard Delaire, 2018
%  
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute GLM';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = '[wip]GLM';
    sProcess.Index       = 1401;
    sProcess.isSeparator = 1;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Movement_correction';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','raw'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    sProcess.options.basis.Comment = 'Basis function';
    sProcess.options.basis.Type    = 'combobox';
    sProcess.options.basis.Value   = {1, {'Hrf','Gamma','BoxCar'}};
    
    sProcess.options.trend.Comment = 'Add a constant trend function ';
    sProcess.options.trend.Type    = 'checkbox';
    sProcess.options.trend.Value   =  1;
    
    sProcess.options.save.Comment = 'Savethe design matrix';
    sProcess.options.save.Type    = 'checkbox';
    sProcess.options.save.Value   =  1;
    
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    
    sProcess.options.fitting.Comment = 'Fitting method';
    sProcess.options.fitting.Type    = 'combobox';
    sProcess.options.fitting.Value   = {1, {'OLS'}};

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
    fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
    if( fitting_choice == 1 )
        Comment=[ Comment ' OLS fit'];
    end    
end

function OutputFiles = Run(sProcess, sInput)

    OutputFiles={};
    
    DataMat = in_bst_data(sInput.FileName);
    basis_choice=sProcess.options.basis.Value{2}(sProcess.options.basis.Value{1});

    % Create the design matrix X : 

    [X,names]=getX(DataMat.Time,DataMat.Events,basis_choice);
    
    
    % Add trend function
    if sProcess.options.trend.Value == 1 
        [C,name]=getTrend(DataMat.Time, 'Constant');
        X=[X C];
        names=[names name];
    end
    
    
    % Check the rank of the matrix
    n_regressor=size(X,2);
    if  rank(X) < n_regressor
        bst_report('Warning', sProcess, sInputs, 'The design matrix is not full-ranked');
    end
    
    % Check the collinearity of the matrix
    
    if cond(X) > 300 
        bst_report('Warning', sProcess, sInputs, [ 'The design matrix is high-correlated : Cond(x)=' num2str(cond(X))] );
    end    

    
    iStudy = sInput.iStudy;

    if sProcess.options.save.Value == 1 
    % Save the results 
        Out_DataMat = db_template('matrixmat');
        Out_DataMat.Value           = X';
        Out_DataMat.Comment     = 'Design Matrix';
        Out_DataMat.Description = names'; % Names of the regressors 
        Out_DataMat.ChannelFlag =  ones(n_regressor,1);   % List of good/bad channels (1=good, -1=bad)
        Out_DataMat.Time        =  DataMat.Time;
        Out_DataMat = bst_history('add', Out_DataMat, 'Design', FormatComment(sProcess));
        OutputFiles{4} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'design_matrix');
        % Save on disk
        save(OutputFiles{4}, '-struct', 'Out_DataMat');
        % Register in database
        db_add_data(iStudy, OutputFiles{4}, Out_DataMat);
    end 
    
    % Solving  B such as Y = XB +e 
    
    %Y= DataMat.F' ;
    B_r=randn(3,60);
    Y=X*B_r;
    Y= Y + + 0.5*randn(size(Y));
    fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
    if( fitting_choice == 1 ) % Use OLS : : \( B= ( X^{T}X)^{-1} X^{T} Y \)  
         B=ols_fit( Y, X );
    else
        bst_report('Error', sProcess, sInputs, 'This method is not implemented yet' );
    end
    
        
    x_hat=X*B;
    
    residual=Y - x_hat ;     
    covE=cov(residual);
    
    for i=1:size(residual,2)     
        covB(:,:,i)=covE(i,i) * inv(transpose(X)*X);
    end
    dfe = size(Y,1) - rank(X);

    
    % Saving the B Matrix.
    
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value           = B';
    Out_DataMat.Comment     = 'B Matrix';
    Out_DataMat.Description = names; % Names of the regressors 
    Out_DataMat.ChannelFlag =  ones(size(B,2),1);   % List of good/bad channels (1=good, -1=bad)
    Out_DataMat = bst_history('add', Out_DataMat, 'GLM computation', FormatComment(sProcess));
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'B_matrix');
    save(OutputFiles{1}, '-struct', 'Out_DataMat');
    db_add_data(iStudy, OutputFiles{1}, Out_DataMat);
    
    % Saving the residual matrix.
    
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value           = residual';
    Out_DataMat.Comment     = 'Residual Matrix';
    %Out_DataMat.Device      = DataMat.Device;
    Out_DataMat.Time        =  DataMat.Time;
    Out_DataMat.Events      =  DataMat.Events;
    Out_DataMat.ChannelFlag =  DataMat.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
    Out_DataMat.DisplayUnits = DataMat.DisplayUnits; 
    
    Out_DataMat = bst_history('add', Out_DataMat, 'GLM computation', FormatComment(sProcess));
    OutputFiles{2} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'residual_matrix');
    save(OutputFiles{2}, '-struct', 'Out_DataMat');
    db_add_data(iStudy, OutputFiles{2}, Out_DataMat);    
    
    % Saving the prediction matrix for debug purpose
    
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value           = x_hat';
    Out_DataMat.Comment     = 'Predict Matrix';
    %Out_DataMat.Device      = DataMat.Device;
    Out_DataMat.Time        =  DataMat.Time;
    Out_DataMat.Events      =  DataMat.Events;
    Out_DataMat.ChannelFlag =  DataMat.ChannelFlag;% List of good/bad channels (1=good, -1=bad)
    Out_DataMat.DisplayUnits = DataMat.DisplayUnits; 
    
    Out_DataMat = bst_history('add', Out_DataMat, 'GLM computation', FormatComment(sProcess));
    OutputFiles{3} = bst_process('GetNewFilename', fileparts(sInput.FileName), ' _matrix');
    save(OutputFiles{3}, '-struct', 'Out_DataMat');
    db_add_data(iStudy, OutputFiles{3}, Out_DataMat);    

end


function B=ols_fit(y,X)

    %B= (transpose(X)*X)\transpose(X)*y;
    %B=X\y;
    
    B= pinv(X'*X)* X'*y;
 
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
        X(:,i) = filter(basis_function(time-time(1)), 1, X(:,i));
        %X(:,i)= conv(X(:,i),basis_function,'same');
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
    if nargin <  2, peakTime = 6; end

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
    if nargin <  2, peakTime    = 6;end
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




