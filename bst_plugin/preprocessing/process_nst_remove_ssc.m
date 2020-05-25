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
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1308;
sProcess.isSeparator = 1;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Optode-separations';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','data'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.nOutputs    = 2;

sProcess.options.separation_threshold_cm.Comment = 'Separation threshold';
sProcess.options.separation_threshold_cm.Type    = 'value';
sProcess.options.separation_threshold_cm.Value   = {1.5, 'cm', 2}; 

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function [OutputFiles_SSC] = Run(sProcess, sInput)

OutputFiles_SSC = {};

separation_threshold_m = sProcess.options.separation_threshold_cm.Value{1} / 100;

    
% Load recordings
if strcmp(sInput.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInput.FileName);
elseif strcmp(sInput.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInput.FileName, [], 1, 1, 'no');
end
    
ChannelMat = in_bst_channel(sInput.ChannelFile);
[nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');

Y= sDataIn.F(nirs_ichans,:)';

model= nst_glm_initialize_model(sDataIn.Time);
model=nst_glm_add_regressors(model,"channel",sInput,'distance', separation_threshold_m);

[B,proj_X] = nst_glm_fit_B(model,Y, 'SVD');
Y= Y - model.X*B;


sDataOut                    = db_template('data');
sDataOut.F                  = sDataIn.F;
sDataOut.F(nirs_ichans,:)   = Y';
sDataOut.Comment            = [sInput.Comment ' SSC'];
sDataOut.ChannelFlag        = sDataIn.ChannelFlag; 
sDataOut.Time               = sDataIn.Time;
sDataOut.DataType           = 'recordings'; 
sDataOut.nAvg               = 1;
if ~isempty(sDataIn.Std)
    sDataOut.Std = sDataIn.Std;
else
    sDataOut.Std = [];
end
sDataOut.ColormapType = [];
sDataOut.Events = sDataIn.Events;
sDataOut.DisplayUnits = sDataIn.DisplayUnits;

% Generate a new file name in the same folder
sStudy = bst_get('Study', sInput.iStudy);
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_lsc');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(sInput.iStudy, OutputFile, sDataOut);
OutputFiles_SSC{1} = OutputFile;

end    