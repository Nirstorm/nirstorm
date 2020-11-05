function varargout = process_nst_motion_correction( varargin )

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
% Authors: Thomas Vincent (2015-2018)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
%TOCHECK: how do we limit the input file types (only NIRS data)?
sProcess.Comment     = 'Motion correction';
sProcess.FileTag     = 'mvt corr';
sProcess.Category    = 'File';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Movement_correction';
sProcess.isSeparator = 0; % add a horizontal bar after the process in
%                             the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options
sProcess.options.option_event_name.Comment = 'Movement event name: ';
sProcess.options.option_event_name.Type    = 'text';
sProcess.options.option_event_name.Value   = '';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFiles = {};

if ~license('test', 'Curve_Fitting_Toolbox')
    bst_error('Curve Fitting Toolbox not available');
    return 
elseif isempty(which('csaps'))
    bst_error(['Curve Fitting Toolbox OK but function csaps not found.<BR>' ...
               'Try refreshing matlab cache using command: rehash toolboxcache']);
    return
end

% Get selected events
event_name =  strtrim(sProcess.options.option_event_name.Value);

for iInput=1:length(sInputs)
    % Load recordings
    if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iInput).FileName);
        events = sDataIn.Events;
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iInput).FileName, 'F');
        events = sDataRaw.F.events;
    end
    
    event = [];
    ievt_mvt = [];
    for ievt=1:length(events)
        if strcmp(events(ievt).label, event_name)
            event = events(ievt);
            ievt_mvt = ievt;
            break;
        end
    end
    if isempty(event)
        warning(['Event "' event_name '" does not exist in file.']);
    end

    if isempty(event) || isempty(event.times) % no marked event
        data_corr = sDataIn.F';
    else
        % Process only NIRS channels
        channels = in_bst_channel(sInputs(iInput).ChannelFile);
        nirs_ichans = channel_find(channels.Channel, 'NIRS');
        data_nirs = sDataIn.F(nirs_ichans, :)';
        prev_negs = any(data_nirs <= 0, 1);
        % Keep track of channels that contained neg values 
        samples = time_to_sample_idx(event.times, sDataIn.Time);
        data_corr = Compute(data_nirs, sDataIn.Time', samples');
        % Fix negative values created by moco
        new_negs = any(data_corr <= 0, 1) & ~prev_negs;
        if any(new_negs)
            warning('Motion correction introduced negative values. Will be fixed by global offset');
            % TODO: use pair-specific offset see /home/tom/Projects/Research/Software/mfipcode/sandbox/nirstorm/mfip_correct_negChan.m
            offset = 2*abs(min(data_corr(:)));
            data_corr = data_corr + offset;
        end
        data_corr_full = sDataIn.F';
        data_corr_full(:, nirs_ichans) = data_corr;
        data_corr = data_corr_full;
    end
    if 0
        nirs_data_full = in_bst(sInputs(iInput).FileName, [], 1, 0, 'no');
        channels = in_bst_channel(sInputs(iInput).ChannelFile);
        [nirs_ichans, tmp] = channel_find(channels.Channel, 'NIRS');
        nirs_chan_flags = zeros(length(channels.Channel), 1);
        nirs_chan_flags(nirs_ichans) = 1;
        [new_ChannelFlag, bad_chan_names] = Compute(nirs_data_full.F', channels, ...
            nirs_data_full.ChannelFlag, ...
            do_remove_neg_channels, max_sat_prop, ...
            invalidate_paired_channels);
        
        % Add bad channels
        % OutputFiles = {sInputs(iInput).FileName};
    end
    
    
    
    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = data_corr';
    sDataOut.Comment      = 'Motion-corrected';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings';
    sDataOut.History      = sDataIn.History;
    sDataOut = bst_history('add', sDataOut, 'process', sProcess.Comment);
    
    sDataOut.nAvg         = 1;
    if ~isempty(ievt_mvt)
        sDataOut.Events       = events([1:(ievt_mvt-1) (ievt_mvt+1):length(events)]);
    else
        sDataOut.Events       = events;
    end
    sDataOut.DisplayUnits = sDataIn.DisplayUnits;
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_motion_corr');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    OutputFiles{iInput} = OutputFile;
end
end


%% ===== Compute =====
function [nirs_sig_corr] = Compute(nirs_sig, t, mvtWin_samples) %#ok<DEFNU>
%% Update the given channel flags to indicate which pairs are to be removed:
%% - negative values
%% - saturating
% TODO: add cirterion on pair distance?
%
% Args
%    - nirs_sig: matrix of double, size: time x nb_channels
%        nirs signals to be corrected
%    - t: array of doublee, size: time x 1
%        time axis of the given nirs signal
%    - mvtWin_samples: matrix of double, size: nb_mvt_events x 2
%        windows of marked movements (starting and ending sample index)
%
% Output:
%    - nirs_sig_corr: matrix of double, size: time x nb_channels
%        Motion-corrected nirs signals
%
%TODO: check that given events are not overlapping and are in
%      ascending onset order.
%TODO: correct negative channels instead. Currently if mvt correction
% produces channels with negative values, they will be filtered out before
% MBLL
%

param.spline.smooth=0.99;

nSamples = size(nirs_sig, 1);
nChan = size(nirs_sig, 2);
nMvt = size(mvtWin_samples, 1);
d = nirs_sig;
%====================================================================
% define normal and mvt segments
%======================================================================


for iM=1:nMvt
    if iM==1
        normalWin_samples(iM,:)=[1 mvtWin_samples(iM,1)-1];
    else
        normalWin_samples(iM,:)=[mvtWin_samples(iM-1,2)+1 mvtWin_samples(iM,1)-1];
    end
    
    if iM==nMvt
        normalWin_samples(iM+1,:)=[mvtWin_samples(iM,2)+1 nSamples];
    end
end
clear iM
nNormalSeg=size(normalWin_samples,1);

%==========================================================================
% Preprocessing MC correction: all channels
%==========================================================================
% F Scholkmann; Physiol.Meas.31 2010

% =============================================================
% segment data
% =============================================================
for iM=1:nMvt
    MA_seg{iM}= d(mvtWin_samples(iM,1):mvtWin_samples(iM,2),:);
end

for iN=1:nNormalSeg
    normal_seg{iN}= d(normalWin_samples(iN,1):normalWin_samples(iN,2),:);
    normal_seg_mean(iN,:)=mean(normal_seg{iN},1);
end
normal_seg_mean_global=mean(normal_seg_mean);

% =============================================================
% estimate the mvt artefact: spline with smooth parameter
% =============================================================
for iM=1:nMvt
    for iChan=1:nChan
        interp_seg{iM}(:,iChan) = fnval(...
            csaps(...
            t(mvtWin_samples(iM,1):mvtWin_samples(iM,2)),...
            MA_seg{iM}(:,iChan),...
            param.spline.smooth),...
            t(mvtWin_samples(iM,1):mvtWin_samples(iM,2)));
    end
end
clear iChan


% =============================================================
% substract the estimated spline from the MA periods
% =============================================================
for iM=1:nMvt
    denoised_seg{iM} = MA_seg{iM} - interp_seg{iM};
    % off=repmat(mean(d([mvtWin_samples(iM,1)-1 mvtWin_samples(iM,2)+1],:)),size(MA_seg{iM},1),1);
end

% =============================================================
% shift all normal segments and replace them in d:
% the method here differ a little from the one proposed by  F Scholkmann
% =============================================================

% -------------------------------------------------------------------------
% % Solution 1: fit to the last sample of the previous segment
% -------------------------------------------------------------------------
%             for iN=2:nNormalSeg
%                 nS_iSeg=size(normal_seg{iN},1);
%                 d(normalWin_samples(iN,1):normalWin_samples(iN,2),:)=normal_seg{iN}-...
%                                                                      repmat(normal_seg_mean(iN,:),nS_iSeg,1)+...
%                                                                      repmat(d(normalWin_samples(iN-1,2),:),nS_iSeg,1);
%             end

% -------------------------------------------------------------------------
% % Solution 2: weighted shift of segments according to the
% number of samples before and after the transition
% -------------------------------------------------------------------------
for iN=1:nNormalSeg-1
    
    nS_iSeg_1=size([1:normalWin_samples(iN,2)]',1);
    nS_iSeg_2=size(d(normalWin_samples(iN,2)+1:end,:),1);
    
    nS_iSeg_1_2=nS_iSeg_1+nS_iSeg_2;
    weight_iSeg_1=(nS_iSeg_1/nS_iSeg_1_2)/0.5;
    weight_iSeg_2=(nS_iSeg_2/nS_iSeg_1_2)/0.5;
    
    val_seg_1=d(normalWin_samples(iN,2),:);
    val_seg_2=d(normalWin_samples(iN+1,1),:);
    shift=( val_seg_2-val_seg_1)/2;
    
    d(1:normalWin_samples(iN,2),:)=d(1:normalWin_samples(iN,2),:)+...
        repmat(weight_iSeg_2*shift,normalWin_samples(iN,2),1);
    
    d(normalWin_samples(iN,2)+1:end,:)=  d(normalWin_samples(iN,2)+1:end,:) -...
        repmat(weight_iSeg_1*shift,size(d(normalWin_samples(iN,2)+1:end,:),1),1);
end



% =============================================================
% add a  constant offset  into the MA residuals to fit previous
% and next normal segments ( interp function is used so if we
% choose another method to shift segments we can add linear
% offset instead)
% =============================================================
for iM=1:nMvt
    off=[];
    for iC=1:nChan
        off(:,iC)=interp1([mvtWin_samples(iM,1)-1;mvtWin_samples(iM,2)+1],...
            [d(mvtWin_samples(iM,1)-1,iC);d(mvtWin_samples(iM,2)+1,iC)],...
            [mvtWin_samples(iM,1)-1:mvtWin_samples(iM,2)+1]',...
            'linear');
    end
    denoised_seg{iM} = denoised_seg{iM} + off(2:end-1,:);
end

% =============================================================
% replace MA segments in d
% =============================================================
for iM=1:nMvt
    d(mvtWin_samples(iM,1):mvtWin_samples(iM,2),:)=denoised_seg{iM};
end

% =============================================================
% offset <0 channels
% =============================================================
% [d ~]=mfip_correct_negChan(d,ml);

nirs_sig_corr = d;
end

function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end