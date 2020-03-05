function varargout = process_nst_extract_ssc( varargin )

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
% Authors: Thomas Vincent (2019)
%
%
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Split short separation channels';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1001;
sProcess.isSeparator = 0;
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

function [OutputFiles_SSC] = Run(sProcess, sInputs)

OutputFiles_SSC = {};

separation_threshold_m = sProcess.options.separation_threshold_cm.Value{1} / 100;

for iInput=1:length(sInputs)
    
    % Load recordings
    if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iInput).FileName);
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
    end
    
    ChannelMat = in_bst_channel(sInputs(iInput).ChannelFile);
    % nirs_ichans = channel_find(ChannelMat.Channel, 'NIRS');
    % other_chans = true(size(sDataIn.ChannelFlag));
    % other_chans(nirs_ichans) = false;
    
    separations = process_nst_separations('Compute', ChannelMat.Channel);
    if isempty(separations)
        warning(sprintf('Separations could not be computed for %s', sInputs.FileName));
        continue;
    end
    
    lsc_chans = separations > separation_threshold_m;
    ssc_chans = ~lsc_chans & ~isnan(separations);
    
    % Create new condition because channel definition is different from original one
    cond_name = sInputs(iInput).Condition;
    if length(cond_name)>=4 && strcmp(cond_name(1:4), '@raw')
        cond_name = cond_name(5:end);
    end
    iStudyOut = db_add_condition(sInputs(iInput).SubjectName, [cond_name, '_SSC']);
    sStudyOut = bst_get('Study', iStudyOut);
        
    % Save SSC
    selected_chans = ssc_chans; % ssc_chans | other_chans;
    
    ChannelOut = ChannelMat;
    ChannelOut.Channel = ChannelOut.Channel(selected_chans);
    ChannelOut.Comment = sprintf('NIRS short channels (%d)', sum(selected_chans));
        % Save channel definition
    [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudyOut);
    db_set_channel(iChannelStudy, ChannelOut, 0, 0);
    
    sDataOut = db_template('data');
    sDataOut.F            = sDataIn.F(selected_chans, :);
    sDataOut.Comment      = [sInputs(iInput).Comment ' SSC'];
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag(selected_chans); 
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = 1;
    if ~isempty(sDataIn.Std)
        sDataOut.Std = sDataIn.Std(selected_chans, :);
    else
        sDataOut.Std = [];
    end
    sDataOut.ColormapType = [];
    sDataOut.Events = sDataIn.Events;
    sDataOut.DisplayUnits = sDataIn.DisplayUnits;

    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudyOut.FileName), 'data_lsc');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudyOut, OutputFile, sDataOut);
    OutputFiles_SSC{iInput} = OutputFile;
end

end
