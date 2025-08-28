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
function sProcess = GetDescription() 
    % Description the process
    sProcess.Comment     = 'GLM - MFX group t-test';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'GLM'};
    sProcess.Index       = 1605;
    sProcess.isSeparator = 0;
    sProcess.Description = '';
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
function Comment = FormatComment(sProcess)  
    Comment = sProcess.Comment; 
end

function OutputFiles = Run(sProcess, sInputs) 
    OutputFiles={};
    
    output_condition = sProcess.options.output_condition.Value;
    output_subject = sProcess.options.output_subject.Value;

    nb_inputs = length(sInputs);
    sConStat = in_bst_data(sInputs(1).FileName);
    con_name = sConStat.contrast_name;
    
    if isfield(sConStat, 'SurfaceFile')
        surface_data = 1;
        nb_positions = size(sConStat.ImageGridAmp, 1);
        all_cons = zeros(nb_inputs, nb_positions);
        all_cons_var = zeros(nb_inputs, nb_positions);
        for iInput=1:nb_inputs
            in_data = in_bst_data(sInputs(iInput).FileName);
            assert(strcmp(con_name, in_data.contrast_name));
            all_cons(iInput, :) = in_data.ImageGridAmp;
            all_cons_var(iInput, :) = in_data.Std.^2;
        end
        n_chan=size(all_cons,2);
        mask=find(~all(all_cons==0,1)); % Only keep channel in the field of view
        all_cons=all_cons(:,mask);
        all_cons_var=all_cons_var(:,mask);
            
    else
        
        % Warnings : This does not check that the montage is consistent
        % across subject 
        surface_data = 0;
        [ChannelMat,ChannelList]=getCommonChannels(sInputs);
        
        n_chan=length(ChannelList);
        all_cons = zeros(nb_inputs, n_chan);
        all_cons_var = zeros(nb_inputs, n_chan);
        for iInput=1:nb_inputs
            in_data=getDataForChannels(sInputs(iInput),ChannelList);
            assert(strcmp(con_name, in_data.contrast_name));
            all_cons(iInput, :) = in_data.F;
            all_cons_var(iInput, :) = in_data.Std.^2;
        end

    end
    
    con_avg = mean(all_cons, 1);
    gp_residuals = all_cons - repmat(con_avg, nb_inputs, 1);
    edf = nb_inputs - 1;
    inter_subj_var = sum(gp_residuals.^2) ./ edf;
    gp_var = all_cons_var + repmat(inter_subj_var, nb_inputs, 1);
    t_stat_value = sum(all_cons ./ gp_var) ./ sqrt(sum(1./gp_var));
    t_stat_value(isnan(t_stat_value)) = 0;
    
    % Calculate p-values from t-values
    p_value = process_test_parametric2('ComputePvalues', t_stat_value, edf, 't', ...
                                 sProcess.options.tail.Value );
    
    p=zeros(1,n_chan);
    t_stat=zeros(1,n_chan);
    df_map = zeros(1,n_chan);
    
    if surface_data
        p(mask)=p_value;
        t_stat(mask)=t_stat_value;
        df_map(mask)=edf;
        
    else
        p(:)=p_value;
        t_stat(:)=t_stat_value;
        df_map(sum(gp_var, 1) > 0) = edf;    
    end    
                             
    % === OUTPUT STRUCTURE ===
    % Initialize output structure

    sOutput = db_template('statmat');
    sOutput.pmap         = p';
    sOutput.tmap         = t_stat';
    sOutput.contrast_name = con_name;
    sOutput.df = df_map';
    
    if surface_data 
        sOutput.HeadModelType = 'surface';
        sOutput.SurfaceFile = sConStat.SurfaceFile;
        sOutput.ChannelFlag = [];
        sOutput.Type = 'results';

    else
        sOutput.ChannelFlag = ones(1,n_chan);
        sOutput.Options.SensorTypes = 'NIRS'; 
        sOutput.Type = 'data';

    end
    
    sOutput.Correction   = 'no';    
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
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
    
    if surface_data
        OutputFiles{1} = bst_process('GetNewFilename', fileparts(sStudy.FileName), 'presults_group_ttest');
        bst_save(OutputFiles{1}, sOutput, 'v6');
        db_add_data(iStudy, OutputFiles{1}, sOutput);
    else 
         % Save channel definition
        [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
        db_set_channel(iChannelStudy, ChannelMat, 0, 0);   
        OutputFiles{1} = bst_process('GetNewFilename', fileparts(sStudy.FileName), 'pdata_group_ttest');
        bst_save(OutputFiles{1}, sOutput, 'v6');
        db_add_data(iStudy, OutputFiles{1}, sOutput);
        
    end


end

function [ChannelMat,ChannelList]=getCommonChannels(sInputs)
    
    nb_inputs=length(sInputs);
    
    % Get the list of the name of commom channel 
    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    ChannelList=getChannelsNames(ChannelMat);
    for iInput=2:nb_inputs
        ChannelMat = in_bst_channel(sInputs(iInput).ChannelFile);
        tmpChannelList=getChannelsNames(ChannelMat);
        % Do not modify the order of the channels
        ChannelList=intersect(ChannelList,tmpChannelList,'stable');
   
    end
    
    
    % Get the list of the commom channel 
    n_channel=length(ChannelList);
    i_chan=1;
    for i_old_chan=1:length(ChannelMat.Channel)
        if ~isempty( intersect( ChannelMat.Channel(i_old_chan).Name,ChannelList)) 
            new_channel(i_chan)=ChannelMat.Channel(i_old_chan);
            i_chan=i_chan+1;
        end
    end
    
    % Update the channel matrix
    ChannelMat.Channel=new_channel;
    ChannelMat.Comment=['NIRS-BRS sensors (' int2str(n_channel) ')'];
end


function ChannelsNames=getChannelsNames(channelMat)
    [nirs_ichans, tmp] = channel_find(channelMat.Channel, 'NIRS');
    ChannelsNames={};
    for i_chan=nirs_ichans
        ChannelsNames{end+1}=channelMat.Channel(i_chan).Name;
    end    
end

function data=getDataForChannels(sInput,ChannelList)
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    channel_name=getChannelsNames(ChannelMat);
    
    data = in_bst_data(sInput.FileName);
    n_chan=length(ChannelList);
    
    F=zeros(n_chan,1);
    contrast_std=zeros(1,n_chan);
            
    for i_chan=1:length(channel_name)
        [C,ia,ib]=intersect(channel_name(i_chan),ChannelList);
        if ~isempty(C) 
            F(ib)=data.F(i_chan);
            contrast_std(ib)=data.Std(i_chan);
        end
    end
    data.F=F;
    data.Std=contrast_std;
    data.ChannelFlag=ones(n_chan,1);
end