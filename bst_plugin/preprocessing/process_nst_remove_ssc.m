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
% Authors: Delaire Edouard (2020)
%
%
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Remove superficial noise';
sProcess.Category    = 'File';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1308;
sProcess.isSeparator = 0;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Optode-separations';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data'};

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
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function [OutputFiles_SSC] = Run(sProcess, sInput)

OutputFiles_SSC = {};


    
% Load recordings
if strcmp(sInput.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInput.FileName);
elseif strcmp(sInput.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInput.FileName, [], 1, 1, 'no');
end
    
ChannelMat = in_bst_channel(sInput.ChannelFile);
[nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');


types=unique({ChannelMat.Channel(nirs_ichans).Group});

% Get time window
if isfield(sProcess.options, 'baseline') && isfield(sProcess.options.baseline, 'Value') && iscell(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value) && ~isempty(sProcess.options.baseline.Value{1})
    Time = sProcess.options.baseline.Value{1};
    iBaseline = panel_time('GetTimeIndices', sDataIn.Time, Time);
else
    Time = [];
    iBaseline = 1:length(sDataIn.Time);
end


F= sDataIn.F;
F_resudial= sDataIn.F;

for itype = 1 :length(types)
    nirs_ichans = strcmp( {ChannelMat.Channel.Type},'NIRS') & strcmp( {ChannelMat.Channel.Group},types{itype});
    model       = nst_glm_initialize_model(sDataIn.Time);
    
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
        Y= sDataIn.F(nirs_ichans,:)';
        F(nirs_ichans,:) = Y';
        continue;
    else
        bst_report('Info',    sProcess, sInput, message)
    end    
    model = nst_glm_add_regressors(model, 'constant');
        

    model_fit = model;
    model_fit.X = model_fit.X(iBaseline,:);
    model_fit.time = model_fit.time(iBaseline);



    Y_baseline = sDataIn.F(nirs_ichans,iBaseline)';
    [B,proj_X] = nst_glm_fit_B(model_fit,Y_baseline, 'SVD');
    

    Y = sDataIn.F(nirs_ichans,:)';
    Y = Y - model.X*B;
    F(nirs_ichans,:) = Y';
    F_resudial(nirs_ichans,:) = (model.X*B)';
end    



sDataOut                    = sDataIn;
sDataOut.F                  = F;
% Add time window to the comment
strMethod = '| SSC (mean)';
if isempty(Time)
    Comment = [strMethod, ': [All file]'];
else 
    Comment = [strMethod, sprintf(': [%1.3fs,%1.3fs]', Time(1), Time(2))];
end
sDataOut.Comment            = [sInput.Comment Comment];
sDataOut                    = bst_history('add', sDataOut, 'process_nst_remove_ssc','Remove superficial noise'); 

% Generate a new file name in the same folder
sStudy = bst_get('Study', sInput.iStudy);
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_lsc');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(sInput.iStudy, OutputFile, sDataOut);
OutputFiles_SSC{1} = OutputFile;



sDataOut                    = sDataIn;
sDataOut.F                  = F_resudial;
% Add time window to the comment
strMethod = '| SSC(residual)';
if isempty(Time)
    Comment = [strMethod, ': [All file]'];
else 
    Comment = [strMethod, sprintf(': [%1.3fs,%1.3fs]', Time(1), Time(2))];
end
sDataOut.Comment            = [sInput.Comment Comment];
sDataOut                    = bst_history('add', sDataOut, 'process_nst_remove_ssc','Remove superficial noise'); 

% % Generate a new file name in the same folder
% sStudy = bst_get('Study', sInput.iStudy);
% OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_lsc2');
% sDataOut.FileName = file_short(OutputFile);
% bst_save(OutputFile, sDataOut, 'v7');
% % Register in database
% db_add_data(sInput.iStudy, OutputFile, sDataOut);
% OutputFiles_SSC{2} = OutputFile;


end    
