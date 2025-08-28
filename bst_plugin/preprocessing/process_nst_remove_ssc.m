function varargout = process_nst_remove_ssc( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Delaire Edouard (2020-2023)
%
%
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
% Description the process
sProcess.Comment     = 'Remove superficial noise';
sProcess.Category    = 'Filter';
sProcess.FileTag     = @GetFileTag;
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1308;
sProcess.isSeparator = 0;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/NIRSTORM#Regressing_out_superficial_noise';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','raw'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.nOutputs    = 1;

sProcess.options.SS_chan.Type       = 'radio_linelabel';
sProcess.options.SS_chan.Comment    = {'Based on Names', 'Based on Source-Detector distances','Short-separation channels: '; 'name', 'distance',''};
sProcess.options.SS_chan.Controller = struct('name','name','distance','distance');
sProcess.options.SS_chan.Value      = 'name';

sProcess.options.SS_chan_name.Comment = str_pad('Superfical Channel [coma-separated list]',44);
sProcess.options.SS_chan_name.Type    = 'text';
sProcess.options.SS_chan_name.Value   = '';     
sProcess.options.SS_chan_name.Class   = 'name';

sProcess.options.separation_threshold_cm.Comment = str_pad('Separation threshold',44);
sProcess.options.separation_threshold_cm.Type    = 'value';
sProcess.options.separation_threshold_cm.Value   = {1.5, 'cm', 2}; 
sProcess.options.separation_threshold_cm.Class   = 'distance';

% === Baseline time window
sProcess.options.baseline.Comment = 'Baseline:';
sProcess.options.baseline.Type    = 'baseline';
sProcess.options.baseline.Value   = [];

end

function s = str_pad(s,padsize)
    if (length(s) < padsize)
        s = [repmat('&nbsp;', 1, padsize - length(s)), s];
    end
    s = ['<FONT FACE="monospace">' s '</FONT>'];
end

%% ===== FORMAT COMMENT =====
function [Comment, fileTag] = FormatComment(sProcess)
    Comment = 'SSC (mean)';
    fileTag = 'scr';
end

%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
    [Comment, fileTag] = FormatComment(sProcess);
end

function sInput = Run(sProcess, sInput)


ChannelMat = in_bst_channel(sInput.ChannelFile);
[nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');


types=unique({ChannelMat.Channel(nirs_ichans).Group});

% Get time window
if isfield(sProcess.options, 'baseline') && isfield(sProcess.options.baseline, 'Value') && iscell(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value{1})
    Time = sProcess.options.baseline.Value{1};
    iBaseline = panel_time('GetTimeIndices', sInput.TimeVector, Time);
else
    Time = [];
    iBaseline = 1:length(sInput.TimeVector);
end


F= sInput.A;
F_resudial= sInput.A;

for itype = 1 :length(types)
    nirs_ichans = strcmp( {ChannelMat.Channel.Type},'NIRS') & strcmp( {ChannelMat.Channel.Group},types{itype});
    model       = nst_glm_initialize_model(sInput.TimeVector);
    
    % Include short-seperation channel
    if strcmp(sProcess.options.SS_chan.Value,'distance') % based on distance

        separation_threshold_m = sProcess.options.separation_threshold_cm.Value{1} / 100; %convert to meter
        [model,code,message]=nst_glm_add_regressors(model,'channel',sInput,'distance', separation_threshold_m,types(itype), []);
    elseif strcmp(sProcess.options.SS_chan.Value,'name') % based on name 

        if ~isempty(sProcess.options.SS_chan_name.Value)
            SS_name=strtrim(split(sProcess.options.SS_chan_name.Value,','));
            [model,code,message]=nst_glm_add_regressors(model,'channel',sInput,'name',SS_name',types(itype),[]);
        end    
    end 
    
    if code < 0 
        bst_report('Error',   sProcess, sInput, message)
        Y = sInput.A(nirs_ichans,:)';
        F(nirs_ichans,:) = Y';
        continue;
    else
        bst_report('Info',    sProcess, sInput, message)
    end    
    model = nst_glm_add_regressors(model, 'constant');
        

    model_fit = model;
    model_fit.X = model_fit.X(iBaseline,:);
    model_fit.time = model_fit.time(iBaseline);



    Y_baseline = sInput.A(nirs_ichans,iBaseline)';
    [B,proj_X] = nst_glm_fit_B(model_fit,Y_baseline, 'SVD');
    

    Y = sInput.A(nirs_ichans,:)';
    Y = Y - model.X*B;
    F(nirs_ichans,:) = Y';
    F_resudial(nirs_ichans,:) = (model.X*B)';
end    


sInput.A  = F;
% Add time window to the comment
strMethod = 'SSC (mean)';
if isempty(Time)
    Comment = [strMethod, ': [All file]'];
else 
    Comment = [strMethod, sprintf(': [%1.3fs,%1.3fs]', Time(1), Time(2))];
end

sInput.CommentTag         = Comment;
sInput.HistoryComment     = ['Remove superficial noise : ' strMethod]; 

end    
