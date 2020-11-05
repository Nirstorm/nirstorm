function varargout = process_nst_convert_to_stat( varargin )
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
% Authors: Thomas Vincent 2018
%  
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Convert to stat map';
    sProcess.Category    = 'Stat1';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1902;
    sProcess.isSeparator = 0;
    sProcess.Description = '';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'pdata'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.comment.Comment = 'Text: ';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)  %#ok<DEFNU>
    Comment = sProcess.Comment;
end

function sOutput = Run(sProcess, sInputs) %#ok<DEFNU>
    sOutput = [];
    
    output_comment = sProcess.options.comment.Value;
    if isempty(output_comment)
        error('Empty output comment');
    end
    
    % Load sInputs
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(1).FileName);
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
    end

    % Output as statmap
    sOutput = db_template('statmat');
    sOutput.tmap         = sDataIn.F(:,1);
    sOutput.pmap         = zeros(size(sOutput.tmap)) + 1e-3; %Dummy value 
    sOutput.df           = ones(size(sDataIn.ChannelFlag)) * 10; %Dummy value 
    sOutput.Correction   = 'no';
    sOutput.Type = 'data';
    sOutput.ChannelFlag = sDataIn.ChannelFlag;
    sOutput.Options.SensorTypes = 'NIRS';
    
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    sOutput.nComponents  = 1;
    
    sOutput.Comment = output_comment;
end
