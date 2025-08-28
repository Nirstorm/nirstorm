function varargout = process_nst_separations( varargin )
% process_nst_separation: compute distances between sources and detectors
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
% Authors: Thomas Vincent, 2018
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
% Description the process
sProcess.Comment     = 'Compute separations';
sProcess.Category    = 'File';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1202;
sProcess.isSeparator = 1;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Optode-separations';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','data'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInput)

    OutputFiles = {};

    % Load recordings
    if strcmp(sInput.FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInput.FileName);
    elseif strcmp(sInput.FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInput.FileName, [], 1, 1, 'no');
    end
    
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');

    separations = Compute(ChannelMat.Channel(nirs_ichans)) * 100; %convert to cm
    

    % Get the output condition - create it if it doesn't exist
    if strcmpi(sInput.FileType, 'raw')
      % Create or get the target study
      newCondition = strrep(sInputs.Condition, '@raw', '');
      iStudy = db_add_condition(sInputs.SubjectName, newCondition);
      db_set_channel(iStudy, sInputs.ChannelFile, 2, 0);
    else
      iStudy = sInput.iStudy;
    end
    sStudy = bst_get('Study', iStudy);

    % Save time-series data
    data_out = zeros(size(sDataIn.F, 1), 1);
    data_out(nirs_ichans,:) = separations;
    sDataOut = db_template('data');
    sDataOut.F            = data_out;
    sDataOut.Comment      = 'Separations';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = [1];
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.DisplayUnits = 'cm';
    
    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_chan_dist');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudy, OutputFile, sDataOut);
    OutputFiles{end+1} = OutputFile;
end


function separations = Compute(channels, pair_ids)
% Compute distances between sources and detectors for given channels or pairs.
% Note: output unit is the same as the one of the input.
%
% Inputs:
%    - channels (struct array):
%       Channel brainstorm structure (see db_template('channeldesc'))
%   [- pair_ids ] (2d array of int): size(nb_pairs x 2)
%      List of pairs for which to compute separations.
%      Column 1 contains source ids and column 2 contains detector ids.
%      Must be consistent with given channels.
%
% Outputs:
%    1D array of double containing separations in the same unit as
%    in given channel.
%    If pair_ids is not given, then size(separations) = nb_channels
%    If pair_ids is given, then size(separations) = size(pair_indexes, 1);

if nargin < 2 % pair_indexes not given
        separations = zeros(length(channels), 1);
        for ichan=1:length(channels)
            if strcmp(channels(ichan).Type, 'NIRS')
                separations(ichan) = euc_dist(channels(ichan).Loc(:,1), ...
                                              channels(ichan).Loc(:,2));
            else
                separations(ichan) = nan;
            end
        end
else
    montage_info = nst_montage_info_from_bst_channels(channels);
    separations = zeros(size(pair_ids, 1), 1);
    for ipair=1:size(pair_ids, 1)
        src_idx = montage_info.src_ids == pair_ids(ipair, 1);
        det_idx = montage_info.det_ids == pair_ids(ipair, 2);
        separations(ipair) = euc_dist(montage_info.src_pos(src_idx, :), ...
                                      montage_info.det_pos(det_idx, :));
    end
end
end

function d = euc_dist(p1, p2)
    d = sqrt(sum((p1 - p2).^2));
end
