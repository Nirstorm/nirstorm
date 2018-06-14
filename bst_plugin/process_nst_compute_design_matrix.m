function varargout = process_nst_compute_design_matrix( varargin )
% process_compute_design_matrix: create the disgn matrix related to event
% specified in the file. 
% 
% Further update : adding trend function
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
%    Define regressors based on events convolved by pre-defined basis functions (canonical HRF, HRF + derivatives...)
%    Define regressors based on external measures
%    Orthogonalize regressors of non-interest onto the event-based ones.
%  
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute design matrix';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'GLM';
    sProcess.Index       = 1400;
    sProcess.isSeparator = 1;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Movement_correction';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','raw'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.basis.Comment = 'Basis function';
    sProcess.options.basis.Type    = 'combobox';
    sProcess.options.basis.Value   = {2, {'Hrf','Gamma'}};
    
    sProcess.options.trend.Comment = 'Add a constant trend function ';
    sProcess.options.trend.Type    = 'checkbox';
    sProcess.options.trend.Value   =  1;
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;

    basis_choice=cell2mat(sProcess.options.basis.Value(1));

    if( basis_choice == 1 )
        Comment=[ Comment ',basis=HRF'];
    elseif ( basis_choice == 2)
        Comment=[ Comment ',basis=gamma'];
    else
        Comment=[ Comment ',basis=HRF'];
    end
    if sProcess.options.trend.Value == 1
       Comment=[ Comment ',constant trend'];
    end
end

function OutputFile = Run(sProcess, sInput)
    DataMat = in_bst_data(sInput.FileName);
    
    events=DataMat.Events;
    n_event=length(events);
    n_sample=length(DataMat.Time);
    
    n_regressor=n_event;
    if sProcess.options.trend.Value == 1
        n_regressor=n_regressor+1;
    end
    X=zeros(n_sample,n_regressor);
    
    
    % create the stim matrix    
    for i=1:n_event
       disp(events(i).label)
       n_run=length(events(i).samples);
       for j =  1:n_run
           disp([ events(i).samples(1,j),events(i).samples(2,j)]);
           event=events(i).samples(1,j):events(i).samples(2,j);
           X(event,i)=ones(size(event));
       end
       
    end
    
    basis_choice=cell2mat(sProcess.options.basis.Value(1));
    % 1= HRF, 2=Gamma, else= HRF
    if( basis_choice == 1 )
        basis_function=Canonical( DataMat.Time);
    elseif ( basis_choice == 2)
        basis_function=Gamma(DataMat.Time);
    else
        basis_function=Canonical( DataMat.Time);
    end
    
    for i=1:n_event
        X(:,i) = filter(basis_function, 1, X(:,i));
        %X(:,i)= conv(X(:,i),basis_function,'same');
    end
    
    % Add trend function
    if sProcess.options.trend.Value == 1 
       X(:,n_event+1)= ones(n_sample,1);
    end
    
    % Check the rank of the matrix
    
    if  rank(X) < n_regressor
        bst_report('Warning', sProcess, sInputs, 'The design matrix is not full-ranked');
    end
    
    % Save the results 
    iStudy = sInput.iStudy;
    
    Out_DataMat = db_template('matrixmat');
    Out_DataMat.Value           = X;
 	Out_DataMat.Comment     = 'Design Matrix';
	Out_DataMat.ChannelFlag =  ones(n_regressor,1);   % List of good/bad channels (1=good, -1=bad)
    Out_DataMat.Time        =  DataMat.Time;
    
    Out_DataMat = bst_history('add', Out_DataMat, 'Design', FormatComment(sProcess));

	OutputFile{1} = bst_process('GetNewFilename', fileparts(sInput.FileName), 'design_matrix');
	% Save on disk
	save(OutputFile{1}, '-struct', 'Out_DataMat');
	% Register in database
	db_add_data(iStudy, OutputFile{1}, Out_DataMat);
           
end