function varargout = process_nst_glm_contrast( varargin )
% process_nst_glm_contrast
% Compute a contrast and its variance from a GLM result
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
% Authors: Edouard Delaire, Thomas Vincent 2018
%  
eval(macro_method);
end



%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - 1st level contrast';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'GLM'};
    sProcess.Index       = 1602;
    sProcess.isSeparator = 0;
    sProcess.Description = '';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results'};
    sProcess.OutputTypes = {'data', 'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
   
    sProcess.options.Contrast.Comment = 'Contrast vector';
    sProcess.options.Contrast.Type    = 'text';
    sProcess.options.Contrast.Value   = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment;
    % Comment = [ Comment ' C = ' sProcess.options.Contrast.Value ]; 
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = [];
        
    if isempty(sProcess.options.Contrast.Value)
       error('Empty contrast');
    end
     
    glm_fit = in_bst_data(sInputs(1).FileName);
    if isfield(glm_fit, 'SurfaceFile') && ~isempty(glm_fit.SurfaceFile)
        surface_data = 1;
        B = glm_fit.ImageGridAmp';
    else
        surface_data = 0;
        B = glm_fit.F';
    end
    covB = glm_fit.beta_cov;
    edf = glm_fit.edf;        
    reg_names = glm_fit.reg_names;
    mse_res = glm_fit.mse_residuals;
    
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
            ind = find(strcmp(reg_names,evt_name));
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
    %covB=covB(:,:,3); 
    con_mat = C * B; 
    if size(covB,3) > 1
        nchan=size(covB,3);
        con_cov=zeros(1,nchan);
        con_std=zeros(1,nchan);
        for ichan=1:nchan
            con_cov(ichan) = C * covB(:,:,ichan) * C';
            con_std(ichan) = sqrt(mse_res(ichan) .* con_cov(ichan));
        end    
    else
        con_cov = C * covB * C';
        con_std = sqrt(mse_res .* con_cov);
    end                             
    % Formating a readable comment such as -Rest +Task
    comment =  strrep(sInputs(1).Comment, ' fitted model' , ' | contrast ');
    contrast = '';
    for i=1:nb_regressors
        if ( C(i) < 0)
            if( C(i) == -1 )
                contrast=[  contrast  ' - ' cell2mat(reg_names(i)) ' '];
            else
                contrast=[  contrast num2str(C(i)) ' ' cell2mat(reg_names(i)) ' '];
            end
        elseif ( C(i) > 0 )
            if( C(i) == 1)
                contrast=[  contrast  ' + ' cell2mat(reg_names(i)) ' '];  
            else
                contrast=[  contrast  ' + ' num2str(C(i)) ' ' cell2mat(reg_names(i)) ' '];
            end 
        end     
    end    
    comment = [comment ' '  contrast];
    
%     iStudy = sInputs.iStudy;    
%     [tmp, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
%     [sStudyIntra, iStudyIntra] = bst_get('AnalysisIntraStudy', iSubject);
%     [ChannelFile] = bst_get('ChannelFileForStudy', iStudy);
%     [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudyIntra);
%     db_set_channel(iChannelStudy, ChannelFile, 0, 0);
%     output_fn = bst_process('GetNewFilename', fileparts(sStudyIntra.FileName), 'pdata_ttest_matrix');

    extra_output.contrast_name = contrast;
    extra_output.Std = con_std';
    extra_output.edf = edf;
    extra_output.DisplayUnits = glm_fit.DisplayUnits;
    
    sStudy = bst_get('Study', sInputs(1).iStudy);

    if surface_data
        [sStudy, ResultFile] = nst_bst_add_surf_data(con_mat', 1, [], 'surf_glm_con', comment, ...
                                                     sInputs(1), sStudy, 'GLM - specify contrast', glm_fit.SurfaceFile, ...
                                                     0, extra_output);
        OutputFiles{end+1} = ResultFile;
    else
        %TODO: check channel-space output
        sDataOut = db_template('data');
        sDataOut.F            = con_mat';
        sDataOut.Comment      = comment;
        sDataOut.ChannelFlag  = glm_fit.ChannelFlag;
        sDataOut.Time         = [1];
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = glm_fit.DisplayUnits; 
        sDataOut.History      = glm_fit.History;
        sDataOut              = bst_history('add',  sDataOut,  'compute',  'GLM - specify contrast');
        
        % Add extra fields
        extra_fields = fieldnames(extra_output);
        for ifield = 1:length(extra_fields)
            if ~isfield(sDataOut, extra_fields{ifield})
                sDataOut.(extra_fields{ifield}) = extra_output.(extra_fields{ifield});
            end
        end
        
        % Save to bst database
        glm_fn = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_glm');
        sDataOut.FileName = file_short(glm_fn);
        bst_save(glm_fn, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs(1).iStudy, glm_fn, sDataOut);
        OutputFiles{end+1} = glm_fn;
    end

 end
