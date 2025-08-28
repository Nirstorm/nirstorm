function varargout = process_nst_merge_montage( varargin )

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
    sProcess.Comment     = 'Merge Montage';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1407;
    sProcess.Description = '';
    sProcess.isSeparator = 1; 
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
%     % Definition of the options

    

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 

   OutputFiles = {};
   
   n_source     = 1;
   n_det        = 1;
   new_channels = []; 
   
   sources_pos  = [];
   det_pos      = [];
   
   for iFile = 1:length(sInputs)

        % Load channel file
        mapSources = containers.Map();
        mapDet     = containers.Map();

        
        ChannelMat = in_bst_channel(sInputs(iFile).ChannelFile);
        channels = ChannelMat.Channel;



        for i_chan = 1:length(channels)
            
            if strcmp(channels(i_chan).Type,'NIRS')
                [isrcs, idets, measures, channel_type] = nst_unformat_channel(channels(i_chan).Name);
                if ~mapSources.isKey(num2str(isrcs))
                    mapSources(num2str(isrcs)) = n_source;
                    sources_pos = [sources_pos channels(i_chan).Loc(:,1) ];
                    n_source = n_source +1;
                end   
                if ~mapDet.isKey(num2str(idets))
                    mapDet(num2str(idets)) = n_det;
                    det_pos = [det_pos channels(i_chan).Loc(:,2) ];
                    n_det = n_det +1;
                end   
                tmp = sprintf('S%dD%dWL%d', mapSources(num2str(isrcs)), mapDet(num2str(idets)),measures);
                %disp(sprintf('%s => %s', channels(i_chan).Name, tmp));
                channels(i_chan).Name = tmp;
            end


            if isempty(new_channels)
                new_channels = channels(i_chan);
            else
                new_channels = [new_channels ; channels(i_chan)];
            end    
        end
   end
   
    msg_s_s = detect_conflict(sources_pos, sources_pos, 1, {'Source', 'Source'});
   
    msg_d_d = detect_conflict(det_pos, det_pos , 1 ,{'Detector', 'Detector'});

    msg_s_d = detect_conflict(sources_pos, det_pos, 0,{'Source', 'Detector'});
    
    if ~isempty(msg_s_s)
        bst_report('Warning', sProcess, sInputs, msg_s_s);
    end
    
    if ~isempty(msg_d_d)
        bst_report('Warning', sProcess, sInputs, msg_d_d);
    end
    
    if ~isempty(msg_s_d)
        bst_report('Warning', sProcess, sInputs, msg_s_d);
    end

    [sSubjStudies, ~] = bst_get('StudyWithSubject', sInputs(1).SubjectFile,'intra_subject', 'default_study');
    newCondition = file_unique([sInputs(1).Condition, '_merge'], {sSubjStudies.Name}, 1);
    iStudy = db_add_condition(sInputs(1).SubjectName, newCondition);

    sStudy = bst_get('Study', iStudy);
    
    ChannelMat         = db_template('channel');
    ChannelMat.Channel = new_channels;
    ChannelMat.Comment = sprintf('NIRS-BRS channels (%d)', length(new_channels));
    for iFile = 1:length(sInputs)
        ChannelMat = bst_history('add', ChannelMat, 'merge', ['Merged file: ' file_short(sInputs(iFile).FileName)]);
        tmp =  in_bst_channel(sInputs(iFile).ChannelFile);
        ChannelMat = bst_history('add', ChannelMat, tmp.History,[ 'From' file_short(sInputs(iFile).FileName) ': ' ] );
        if isfield(tmp,'Nirs') && ~isfield(ChannelMat,'Nirs')
            ChannelMat.Nirs = tmp.Nirs;
        end
    end

    % Save channel definition
    [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
    db_set_channel(iChannelStudy, ChannelMat, 2, 0);
    
    separations = process_nst_separations('Compute',ChannelMat.Channel) * 100; %convert to cm
    % Save time-series data
    sDataOut              = db_template('data');
    sDataOut.F            = separations;
    sDataOut.Comment      = 'Separations';
    sDataOut.ChannelFlag  = ones(length(separations),1);
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
    OutputFiles{1} = OutputFile;
    
end

function d = euc_dist(p1, p2)
    d = sqrt(sum((p1 - p2).^2));
end

function msg = detect_conflict(pos1, pos2, isSymetrical, labels)
    msg = '';
    conflict =cell(1, size(pos1,2));
    for i = 1:size(pos1,2)
        D = zeros(1, size(pos2,2));
        if isSymetrical
            n2 = i;
        else
            n2 = size(pos2,2);
        end    
        for j = 1:n2
            D(j) = euc_dist(pos1(:,i), pos2(:,j));
        end
        colapsing = find( D > 0 & D*100 < 1.5);
        if ~isempty(colapsing)
              conflict{i} = union( conflict{i}, colapsing)';  
        end    
    end  

    if any(~cellfun(@isempty,conflict))
        msg = [msg, sprintf('%s - %s conflicts: \n', labels{1},labels{2})];
        for i=1:length(conflict)
            if ~isempty(conflict{i})
                msg = [msg, sprintf('%s %d and %s %s are conflicting \n',labels{1}, i, labels{2}, strjoin(strsplit(num2str(conflict{i})), ' ; '))];
           end 
        end
    end
end   
