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
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1308;
sProcess.isSeparator = 0;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Optode-separations';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','data'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.nOutputs    = 2;

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

F= sDataIn.F;
for itype = 1 :length(types)
    
    model= nst_glm_initialize_model(sDataIn.Time);
    
    % Include short-seperation channel
    if strcmp(sProcess.options.SS_chan.Value,'distance') % based on distance

        separation_threshold_m = sProcess.options.separation_threshold_cm.Value{1} / 100; %convert to meter
        model=nst_glm_add_regressors(model,'channel',sInput,'distance', separation_threshold_m,types(itype));

    elseif strcmp(sProcess.options.SS_chan.Value,'name') % based on name 

        if ~isempty(sProcess.options.SS_chan_name.Value)
            SS_name=split(sProcess.options.SS_chan_name.Value,',');
            model=nst_glm_add_regressors(model,'channel',sInput,'name',SS_name',types(itype));
        end    
    end 
    model = nst_glm_add_regressors(model, "constant");
    
    nirs_ichans = strcmp( {ChannelMat.Channel.Type},'NIRS') & strcmp( {ChannelMat.Channel.Group},types{itype});
    Y= sDataIn.F(nirs_ichans,:)';

    [B,proj_X] = nst_glm_fit_B(model,Y, 'SVD');
    Y= Y - model.X*B;
    F(nirs_ichans,:) = Y';
end    

sDataOut                    = sDataIn;
sDataOut.F                  = F;
sDataOut.Comment            = [sInput.Comment '| SSC'];
sDataOut                    = bst_history('add', sDataOut, 'process_nst_remove_ssc','Remove superficial noise'); 

% Generate a new file name in the same folder
sStudy = bst_get('Study', sInput.iStudy);
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_lsc');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(sInput.iStudy, OutputFile, sDataOut);
OutputFiles_SSC{1} = OutputFile;

end    