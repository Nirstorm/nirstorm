function varargout = process_nst_rename_montage( varargin )

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
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
% Authors: Edouard Delaire, 2026 ; Jean-Eudes Bornert, 2026


eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    sProcess.Comment     = 'Rename Montage';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Optimal Montage'};
    sProcess.Index       = 1103;
    sProcess.Description = '';
    sProcess.isSeparator = 1; 
    
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.montage_panel.Comment = {'panel_nst_rename_montage', 'Edit optodes name:'};
    sProcess.options.montage_panel.Type    = 'editpref';
    sProcess.options.montage_panel.Value   = []; 
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
    OutputFiles = {sInputs.FileName};
    
    panelData = sProcess.options.montage_panel.Value;
    
    if isempty(panelData) || ~isfield(panelData, 'transformations')
        disp('BST> No transformation defined or operation canceled. ');
        return;
    end
    
    user_transformation = panelData.transformations;
    sChannel = in_bst_channel(sInputs.ChannelFile);
    
    for iChannel = 1:length(sChannel.Channel)
        old_channel_name = sChannel.Channel(iChannel).Name;
        
        [isrc, idet, measure, ~] = nst_unformat_channel(old_channel_name);
        source_name = sprintf('S%d', isrc);
        detector_name = sprintf('D%d', idet);
        
        if isfield(user_transformation, source_name)
            source_name = user_transformation.(source_name);
        end
        if isfield(user_transformation, detector_name)
            detector_name = user_transformation.(detector_name);
        end
        
        sChannel.Channel(iChannel).Name = sprintf('%s%sWL%d', source_name, detector_name, measure);
    end
    
    if ~isempty(sChannel.Clusters)
        warning('Removing clusters');
        sChannel.Clusters = [];
    end
    bst_save(file_fullpath(sInputs.ChannelFile), sChannel);
end