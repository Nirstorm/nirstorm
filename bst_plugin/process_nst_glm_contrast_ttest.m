function varargout = process_nst_glm_contrast_ttest( varargin )
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
% Authors: Edouard Delaire, Thomas Vincent 2018
%  
eval(macro_method);
end



%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - contrast t-test';
    sProcess.Category    = 'Stat1';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results'};
    sProcess.OutputTypes = {'pdata', 'presults'}; %TODO sort out output types depending on data type (chan or surf)
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % === tail for the ttest 
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'one+';  
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment;
    % Comment = [ Comment ' C = ' sProcess.options.Contrast.Value ]; 
end

function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    sOutput = [];
    
    con_data = in_bst_data(sInputs(1).FileName);
    if isfield(con_data, 'SurfaceFile')
        surface_data = 1;
        con_mat = con_data.ImageGridAmp';
    else
        con_mat = con_data.F;
    end
    edf = con_data.edf;        
    con_std = con_data.contrast_std;

    t_stat = con_mat ./ con_std ;
    t_stat(con_std==0) = 0;
    p = process_test_parametric2('ComputePvalues', t_stat, edf, 't', ...
                                 sProcess.options.tail.Value );
                             
    comment = 't-stat';
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
    sOutput.pmap         = p';
    sOutput.tmap         = t_stat';
    sOutput.contrast_name = con_data.contrast_name;
    
    df_map = zeros(size(sOutput.tmap));
    df_map(con_std>0) = edf;
    sOutput.df = df_map;
    
    sOutput.Correction   = 'no';
    if surface_data
        sOutput.SurfaceFile = con_data.SurfaceFile;
        sOutput.ChannelFlag = [];
    else
        sOutput.ChannelFlag = ones(1,nb_positions); %TODO: use current channel flags
        sOutput.Options.SensorTypes = 'NIRS';
    end
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    sOutput.nComponents  = 1;
    
    sOutput.Comment = comment;
    
    sOutput = bst_history('add', sOutput, con_data.History, '');
    sOutput = bst_history('add', sOutput, 'ttest computation', comment);
        
end
