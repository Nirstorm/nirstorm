function varargout = process_nst_detrend( varargin )

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
sProcess.Comment     = 'Remove slow fluctuations';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1309;
sProcess.isSeparator = 1;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Optode-separations';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','data'};

sProcess.options.option_period.Comment = 'Miminum Period (0= only linear detrend):';
sProcess.options.option_period.Type    = 'value';
sProcess.options.option_period.Value   = {200, 's', 0};

sProcess.options.option_keep_mean.Comment = 'Keep the mean';
sProcess.options.option_keep_mean.Type    = 'checkbox';
sProcess.options.option_keep_mean.Value   = 0;
    
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.nOutputs    = 1;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInput)

OutputFiles = {};
keep_mean = sProcess.options.option_keep_mean.Value;
    
% Load recordings
if strcmp(sInput.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInput.FileName);
elseif strcmp(sInput.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInput.FileName, [], 1, 1, 'no');
end
    
ChannelMat = in_bst_channel(sInput.ChannelFile);
nirs_ichans = sDataIn.ChannelFlag ~= -1 & strcmpi({ChannelMat.Channel.Type}, 'NIRS')';


Y= sDataIn.F(nirs_ichans,:)';

model= nst_glm_initialize_model(sDataIn.Time);
model=nst_glm_add_regressors(model,"constant");
model=nst_glm_add_regressors(model,"linear");  


if sProcess.options.option_period.Value{1} > 0
   period=sProcess.options.option_period.Value{1};
   model=nst_glm_add_regressors(model,"DCT",[1/sDataIn.Time(end) 1/period],{'LFO'});  
end
       
[B,proj_X] = nst_glm_fit_B(model,Y, 'SVD');

if keep_mean
    Y= Y - model.X(:,2:end)*B(2:end,:);
else    
    Y= Y - model.X*B;
end

%disp(sprintf('Rank X: %d, dim X: %d',rank(model.X),size(model.X,2)))

figure(1); 
for i=1:size(model.X,2)
    subplot(size(model.X,2),1,i)
    plot(model.time,model.X(:,i))
    title(model.reg_names{i})
end 

sDataOut                    = db_template('data');
sDataOut.F                  = sDataIn.F;
sDataOut.F(nirs_ichans,:)   = Y';
sDataOut.Comment            = [sInput.Comment ' detrend'];
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
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_detrend');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(sInput.iStudy, OutputFile, sDataOut);
OutputFiles{1} = OutputFile;

end    