function varargout = process_nst_glm_group_ttest( varargin )
% process_nst_compute_group_ttest
% Compute a t-test using mixed effects
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
% Authors: Thomas Vincent, Edouard Delaire, 2018-2019
%  
eval(macro_method);
end



%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - MFX group t-test';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1403;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-implementation';
    % todo add a new tutorial
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'results'};
    sProcess.OutputTypes = {'pdata', 'presults'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;

%     % === Method for the group analysis
%     sProcess.options.method.Comment  = {'Mean',''; ...
%                                         'mean', ''};
%     sProcess.options.method.Type     = 'radio_linelabel';
%     sProcess.options.method.Value    = 'mean';

    sProcess.options.output_subject.Comment = 'Subject: ';
    sProcess.options.output_subject.Type    = 'text';
    sProcess.options.output_subject.Value   = 'Group analysis';

    sProcess.options.output_condition.Comment = 'Condition: ';
    sProcess.options.output_condition.Type    = 'text';
    sProcess.options.output_condition.Value   = 'GLM';
    
    % === TAIL FOR THE TEST STATISTIC
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'two';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment; 
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles={};
    
    output_condition = sProcess.options.output_condition.Value;
    output_subject = sProcess.options.output_subject.Value;

    nb_inputs = length(sInputs);
    sConStat = in_bst_data(sInputs(1).FileName);
    con_name = sConStat.contrast_name;
    
    if isfield(sConStat, 'SurfaceFile')
        nb_positions = size(sConStat.ImageGridAmp, 1);
        all_cons = zeros(nb_inputs, nb_positions);
        all_cons_var = zeros(nb_inputs, nb_positions);
        for iInput=1:nb_inputs
            in_data = in_bst_data(sInputs(iInput).FileName);
            assert(strcmp(con_name, in_data.contrast_name));
            all_cons(iInput, :) = in_data.ImageGridAmp;
            all_cons_var(iInput, :) = in_data.contrast_std.^2;
        end
    else
        %TODO
        error('TODO');
    end
    
    con_avg = mean(all_cons, 1);
    gp_residuals = all_cons - repmat(con_avg, nb_inputs, 1);
    edf = nb_inputs - 1;
    inter_subj_var = sum(gp_residuals.^2) ./ edf;
    gp_var = all_cons_var + repmat(inter_subj_var, nb_inputs, 1);
    t_stat = sum(all_cons ./ gp_var) ./ sqrt(sum(1./gp_var));
    t_stat(isnan(t_stat)) = 0;
    
    % Calculate p-values from t-values
    p = process_test_parametric2('ComputePvalues', t_stat, edf, 't', ...
                                 sProcess.options.tail.Value );
    
    if ~isfield(sConStat, 'SurfaceFile') || isempty(sConStat.SurfaceFile) 
        %TODO
        error('TODO');
        [ChannelFile] = bst_get('ChannelFileForStudy', iStudy);
        [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudyIntra);
        db_set_channel(iChannelStudy, ChannelFile, 0, 0);
    end
    
    % === OUTPUT STRUCTURE ===
    % Initialize output structure

    sOutput = db_template('statmat');
    sOutput.pmap         = p';
    sOutput.tmap         = t_stat';
    sOutput.contrast_name = con_name;
    
    df_map = zeros(size(sOutput.tmap));
    df_map(sum(gp_var, 1) > 0) = edf;
    sOutput.df = df_map;
    
    if ~isfield(sConStat, 'SurfaceFile') || isempty(sConStat.SurfaceFile)
        sOutput.ChannelFlag = ones(1,n_chan);
        sOutput.Options.SensorTypes = 'NIRS';
    else
        sOutput.HeadModelType = 'surface';
        sOutput.SurfaceFile = sConStat.SurfaceFile;
        sOutput.ChannelFlag = [];
    end
    sOutput.Correction   = 'no';    
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    
    sOutput.Type = 'results';
    sOutput.nComponents = 1; %TODO: check
    
    comment = sprintf('GLM group contrast (n=%d df):%s', nb_inputs, con_name);
    switch sProcess.options.tail.Value 
        case {'one-'}
             comment = [ comment ' < 0 '];      
        case {'two'}
             comment = [ comment ' <> 0 ']; 
        case { 'one+'}
             comment = [ comment ' > 0 ']; 
    end  
    
    sOutput.Comment = comment;
    sOutput = bst_history('add', sOutput, 'Group stat', comment);
    
    [sSubject, iSubject] = bst_get('Subject', output_subject, 1);
    if isempty(sSubject)
        [sSubject, iSubject] = db_add_subject(output_subject, [], 1, 0);
    end
    
    iStudy = db_add_condition(output_subject, output_condition);
    sStudy = bst_get('Study', iStudy);
    
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sStudy.FileName), 'presults_group_ttest');
    bst_save(OutputFiles{1}, sOutput, 'v6');
    db_add_data(iStudy, OutputFiles{1}, sOutput);
end
