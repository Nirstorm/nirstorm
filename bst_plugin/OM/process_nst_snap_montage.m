function varargout = process_nst_snap_montage( varargin )

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
% Authors: Edouard Delaire, 2020; Thomas Vincent, 2015-2019


eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    sProcess.Comment     = 'Snap Montage to cap';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Optimal Montage'};
    sProcess.Index       = 1103;
    sProcess.Description = '';
    sProcess.isSeparator = 0; 
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % Definition of the options
    % Description of the process
    sProcess.options.user_input.Comment = 'Condition with the cap';
    sProcess.options.user_input.Type    = 'text';
    sProcess.options.user_input.Value    = '';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 

    OutputFiles     = {sInputs.FileName};
    cap_condition   = sProcess.options.user_input.Value;
    
    if length(unique({sInputs.SubjectName})) > 1
        bst_error('Cannot process multiple subjects in process snap to cap')
        return
    end
    
    sStudy = bst_get('StudyWithCondition', sprintf('%s/%s', sInputs(1).SubjectName, cap_condition));
    sCap =  in_bst_channel(sStudy.Channel.FileName);
    cap_position = [sCap.Channel.Loc];
    

    unique_input_conditions = unique({sInputs.Condition});
    for iCondition =  1:length(unique_input_conditions)
        iFile = find(strcmp({sInputs.Condition}, unique_input_conditions{iCondition}), 1);

        sChannel = in_bst_channel(sInputs(iFile).ChannelFile);
        for iChannel = 1:length(sChannel.Channel)
            
            src_loc = sChannel.Channel(iChannel).Loc(:, 1);
            [~, idx_source] = min(nst_pdist(src_loc', cap_position'));
            sChannel.Channel(iChannel).Loc(:, 1) = cap_position(:, idx_source);
    
    
            det_loc = sChannel.Channel(iChannel).Loc(:, 2);
            [~, idx_det] = min(nst_pdist(det_loc', cap_position'));
            sChannel.Channel(iChannel).Loc(:, 2) = cap_position(:, idx_det);
    
            sChannel.Channel(iChannel).Comment = sprintf('%s => %s', sCap.Channel(idx_source).Name, sCap.Channel(idx_det).Name);
        end
        
        sChannel = bst_history('Add', sChannel, 'Compute', sprintf('Snap position to cap: %s', cap_condition));
        bst_save(file_fullpath(sInputs(iFile).ChannelFile), sChannel);
    end

end
