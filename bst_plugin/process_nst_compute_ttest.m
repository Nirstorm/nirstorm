function varargout = process_nst_compute_ttest( varargin )
% process_nst_compute_ttest
% Compute a ttest according to a spm-like constrast vector using 
% t = cB / sqrt( c Cov(B) c^T )  
%
% B, cov(B) and the corresponding degrree of freedom are estimated inprocess_nst_compute_glm.
% 
% 
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
    sProcess.Comment     = 'GLM - intra subject contrast';
    sProcess.Category    = 'Stat1'; %'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix', 'matrix'};
    sProcess.OutputTypes = {'data', 'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
   
    sProcess.options.Contrast.Comment = 'Contrast vector';
    sProcess.options.Contrast.Type    = 'text';
    sProcess.options.Contrast.Value   = '';
    
    % === tail for the ttest 
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'one+';

    
    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
    % Comment = [ Comment ' C = ' sProcess.options.Contrast.Value ]; 
end

function sOutput = Run(sProcess, sInputs)
    sOutput = [];
        
    if isempty(sProcess.options.Contrast.Value)
       error('Empty contrast');
    end
    
    glm_fit = in_bst_matrix(sInputs(1).FileName);
    B = glm_fit.beta;
    covB = glm_fit.beta_cov;
    df = glm_fit.df;        
       
    nb_regressors = size(B,1);
    nb_positions = size(B,2);
    
    % exctract the constrast vector. 
    if( strcmp( sProcess.options.Contrast.Value(1),'[') && strcmp( sProcess.options.Contrast.Value(end),']') )
        % The constrast vector is in a SPM-format : 
        % sProcess.options.Contrast.Value = '[1,0,-1]'
        C=strsplit( sProcess.options.Contrast.Value(2:end-1),',');
        C=str2num(cell2mat(C));
    else 
        C=[];
        % Parse the input 
        % Accepted form are : 'X event1 Y event2' 
        % where X,Y is a sign(+,-), a signed number(-1,+1) or a decimal(+0.5,-0.5)
        % event12 are the name of regressor present in the design matrix
        % A valid input can be '- rest + task' ( spaces are important)
        
        expression='[-+]((\d+.\d+)|((\d)+)|)\s+\w+';
        [startIndex,endIndex] = regexp( sProcess.options.Contrast.Value , expression ); 

        for(i=1:length(startIndex))
            word=sProcess.options.Contrast.Value(startIndex(i):endIndex(i)); % can be '-rest','+rest'..
            
            [evt_ind_start,evt_ind_end]=regexp(word, '\s+\w+' );            
            evt_name=word(evt_ind_start+1:evt_ind_end);

            
            % Find the weight of the regressor 
            if strcmp(word(1:evt_ind_start-1),'+')
                evt_coef=1;
            elseif strcmp(word(1:evt_ind_start-1),'-')
                evt_coef=-1;
            else
                evt_coef=str2double(word(1:evt_ind_start));
            end
            
            
            %Find the position of the regressor            
            ind = find(strcmp(glm_fit.Description,evt_name));
            if( isempty(ind) || ~isnumeric(evt_coef) )
               bst_report('Error', sProcess, sInputs, [ 'Event ' evt_name ' has not been found']);
               return;
            end
            
            C(ind)=evt_coef;
        end

        if isempty(C)
            bst_report('Error', sProcess, sInputs, 'The format of the constrast vector (eg [-1 1] ) is not recognized');
            return
        end
   end
    
    % Add zero padding for the trend regressor 
    if length(C) < nb_regressors
       C= [C zeros(1, nb_regressors - length(C)) ]; 
    end    
     
    con_mat = C * B;
    t_stat = zeros(1, nb_positions);
    
    %TODO: remove loop over positions
    for i = 1:nb_positions
        if all(covB(:,:,i)~=0) % skip positions where signal was null TODO: use proper masking
            t_stat(i) = con_mat(i) / sqrt( C*covB(:,:,i)*transpose(C) ) ;
        end
    end
    
    p = process_test_parametric2('ComputePvalues', t_stat, df, 't', ...
                                 sProcess.options.tail.Value );
                             
   % Formating a readable comment such as -Rest +Task
    comment='T-test : ';
    contrast='';
    for i=1:nb_regressors
        if ( C(i) < 0)
            if( C(i) == -1 )
                contrast=[  contrast  ' - ' cell2mat(glm_fit.Description(i)) ' '];
            else
                contrast=[  contrast num2str(C(i)) ' ' cell2mat(glm_fit.Description(i)) ' '];
            end
        elseif ( C(i) > 0 )
            if( C(i) == 1)
                contrast=[  contrast  ' + ' cell2mat(glm_fit.Description(i)) ' '];  
            else
                contrast=[  contrast  ' + ' num2str(C(i)) ' ' cell2mat(glm_fit.Description(i)) ' '];
            end 
        end     
    end    
    comment=[comment contrast];
    switch sProcess.options.tail.Value 
        case {'one-'}
             comment=[ comment ' < 0 '];      
        case {'two'}
             comment=[ comment ' <> 0 ']; 
        case { 'one+'}
             comment=[ comment ' > 0 ']; 
    end    
    

    
%     iStudy = sInputs.iStudy;    
%     [tmp, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
%     [sStudyIntra, iStudyIntra] = bst_get('AnalysisIntraStudy', iSubject);
%     [ChannelFile] = bst_get('ChannelFileForStudy', iStudy);
%     [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudyIntra);
%     db_set_channel(iChannelStudy, ChannelFile, 0, 0);
%     output_fn = bst_process('GetNewFilename', fileparts(sStudyIntra.FileName), 'pdata_ttest_matrix');
    
    % Output of statmap
    sOutput = db_template('statmat');
    sOutput.pmap         = reshape(p, nb_positions, 1); 
    sOutput.tmap         = reshape(t_stat, nb_positions, 1);
    sOutput.df           = ones(nb_positions, 1) * df;
    sOutput.Correction   = 'no';
    if isfield(glm_fit, 'SurfaceFile')
        sOutput.Type = 'results';
        sOutput.SurfaceFile = glm_fit.SurfaceFile;
        sOutput.ChannelFlag = [];
    else
        sOutput.Type = 'data';
        sOutput.ChannelFlag = ones(1,nb_positions); %TODO: use current channel flags
        sOutput.Options.SensorTypes = 'NIRS';
    end
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    sOutput.nComponents  = 1;
    
    sOutput.Comment = comment;
    
    sOutput = bst_history('add', sOutput, glm_fit.History, '');
    sOutput = bst_history('add', sOutput, 'ttest computation', comment);
    
%     OutputFiles{1} = output_fn;
%     save(OutputFiles{1}, '-struct', 'sOutput');
%     db_add_data(iStudyIntra, OutputFiles{1}, sOutput);

%         % Saving the cB Matrix
% TODO: better save as either channel-space map or cortical map
%     sOutput_b = db_template('matrixmat');
%     sOutput_b.F           = con_mat';
%     sOutput_b.Comment     = [contrast ' B' ];
%     sOutput_b.Description = contrast;  
%     sOutput_b.ChannelFlag =  B.ChannelFlag;
%     sOutput_b.Time         = [1];
%     sOutput_b.DataType     = 'recordings';
%     sOutput_b.nAvg         = 1;
%     sOutput_b.DisplayUnits = 'mmol.l-1'; %TODO: check scaling
% 
%     sOutput_b = bst_history('add', sOutput_b, B.History, '');
%     sOutput_b = bst_history('add', sOutput_b, 'Subject Stat ', comment);
%     
%     OutputFiles{2} = bst_process('GetNewFilename', fileparts(sStudyIntra.FileName), 'data_cbeta_matrix');
%     save(OutputFiles{2}, '-struct', 'sOutput_b');
%     db_add_data(iStudyIntra, OutputFiles{2}, sOutput_b);
end
