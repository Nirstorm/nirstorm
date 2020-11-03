function varargout = process_nst_detect_bad( varargin )

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
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    %TOCHECK: how do we limit the input file types (only NIRS data)?
    sProcess.Comment     = 'Detect bad channels';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1301;  
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Bad_channel_tagging';
    sProcess.isSeparator = 0; 
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % Definition of the options
    sProcess.options.text1.Comment   = '<b>Data based detection</b>'; 
    sProcess.options.text1.Type    = 'label';
    
    sProcess.options.option_sci.Comment = 'Channel rejection by Scalp Coupling Index';
    sProcess.options.option_sci.Type    = 'checkbox';
    sProcess.options.option_sci.Value   = 0;
    sProcess.options.option_sci.Controller='sci';
    
    sProcess.options.sci_threshold.Comment = 'SCI threshold:';
    sProcess.options.sci_threshold.Type    = 'value';
    sProcess.options.sci_threshold.Value   = {80, '%', 0};
    sProcess.options.sci_threshold.Class='sci';
    
    sProcess.options.option_coefficient_variation.Comment = 'Channel rejection by coefficient variation';
    sProcess.options.option_coefficient_variation.Type    = 'checkbox';
    sProcess.options.option_coefficient_variation.Value   = 0;
    sProcess.options.option_coefficient_variation.Controller='variation';
    
    sProcess.options.coefficient_variation.Comment = 'Channel rejection by coefficient variation';
    sProcess.options.coefficient_variation.Type    = 'value';
    sProcess.options.coefficient_variation.Value   = {10, '%', 0};
    sProcess.options.coefficient_variation.Class='variation';
    
    
    sProcess.options.option_remove_saturating.Comment = 'Remove saturating channels';
    sProcess.options.option_remove_saturating.Type    = 'checkbox';
    sProcess.options.option_remove_saturating.Value   = 0;
    sProcess.options.option_remove_saturating.Controller='saturation';
    
    sProcess.options.option_max_sat_prop.Comment = 'Maximum proportion of saturating points';
    sProcess.options.option_max_sat_prop.Type    = 'value';
    sProcess.options.option_max_sat_prop.Value   = {10, '%', 0};
    sProcess.options.option_max_sat_prop.Class   = 'saturation';

    sProcess.options.option_min_sat_prop.Comment = 'Maximum proportion of flooring points';
    sProcess.options.option_min_sat_prop.Type    = 'value';
    sProcess.options.option_min_sat_prop.Value   = {10, '%', 0};
    sProcess.options.option_min_sat_prop.Class   = 'saturation';

    sProcess.options.text2.Comment   = '<b>Montage based detection</b>'; 
    sProcess.options.text2.Type    = 'label';
    
    sProcess.options.option_separation_filtering.Comment = 'Filter channels based on separation';
    sProcess.options.option_separation_filtering.Type    = 'checkbox';
    sProcess.options.option_separation_filtering.Value   = 0;
    sProcess.options.option_separation_filtering.Controller   = 'separation';
    
    sProcess.options.option_separation.Comment = 'Acceptable separation: ';
    sProcess.options.option_separation.Type    = 'range';
    sProcess.options.option_separation.Value   = {[0, 5], 'cm', 2};
    sProcess.options.option_separation.Class   = 'separation';
    
    sProcess.options.text3.Comment   = '<b>Other</b>'; 
    sProcess.options.text3.Type    = 'label';
    
    sProcess.options.auxilary_signal.Comment   = 'Auxilary measurment:';
    sProcess.options.auxilary_signal.Type    = 'combobox';
    sProcess.options.auxilary_signal.Value   = {1, {'Keep all','Remove flat','Remove all'}};
    
    sProcess.options.option_keep_unpaired.Comment = 'Keep unpaired channels';
    sProcess.options.option_keep_unpaired.Type    = 'checkbox';
    sProcess.options.option_keep_unpaired.Value   = 0;
    
    %TODO: scalp contact index and outlier detection mentioned by Zhengchen.
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
        
    % Load channel file
    ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
    
    % Load recordings
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sData = in_bst_data(sInputs(1).FileName);
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sData = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
    end
    nb_channel=sum(strcmpi({ChanneMat.Channel.Type}, 'NIRS'));
    [new_ChannelFlag, bad_chan_names,msg] = Compute(sData, ChanneMat, ...
                                                sProcess.options);
    
    for i=1:size(msg,1)
        fraction= round(100*msg{i,2}/nb_channel);
        if fraction > 80  % send an error if more than 80% of the channels are marked bad
            bst_report('Error',    sProcess, sInputs, sprintf('Number of %s: %d (%d%% of the channels)', msg{i,1},msg{i,2},fraction))
        elseif fraction> 20 % send a warning if more than 20% of the channels are marked bad
            bst_report('Warning',    sProcess, sInputs, sprintf('Number of %s: %d (%d%% of the channels)', msg{i,1},msg{i,2},fraction))
        elseif msg{i,2} > 0
            bst_report('Info',    sProcess, sInputs, sprintf('Number of %s: %d', msg{i,1},msg{i,2}))
        else 
            bst_report('Info',    sProcess, sInputs, sprintf('No %s', msg{i,1}))
        end
    end    
                                            
    % Add bad channels
    tree_set_channelflag({sInputs.FileName}, 'AddBad', bad_chan_names);
    OutputFiles = {sInputs.FileName};
end


%% ===== Compute =====
function [channel_flags, removed_channel_names,msg] = Compute(sData, channel_def, options)
%% Update the given channel flags to indicate which pairs are to be removed:
%% - negative values
%% - saturating
%% - too long or too short separation
%
% Args
%    - nirs_sig: matrix of double, size: time x nb_channels 
%        nirs signals to be filtered
%    - channel_def: struct
%        Defintion of channels as given by brainstorm
%        Used fields: Nirs.Wavelengths, Channel
%    - channels_flags: array of int, size: nb_channels
%        channel flags to update (1 is good, 0 is bad)
%   [- do_remove_neg_channels: boolean], default: 1
%        actually remove pair where at least one channel has negative
%        values
%   [- max_sat_prop: double between 0 and 1], default: 1
%        maximum proportion of saturating values.
%        If 1 then all time series can be saturating -> ignore
%        If 0.2 then if more than 20% of values are equal to the max
%        the pair is discarded.
%   [- min_sat_prop: double between 0 and 1], default: 1
%        maximum proportion of flooring values.
%        If 1 then all time series can be flooring (equal to lowest value) -> ignore
%        If 0.2 then if more than 20% of values are equal to the min value
%        the pair is discarded.
%   [- max_separation_cm: positive double], default: 10
%        maximum optode separation in cm.
%   [- min_separation_cm: positive double], default: 0
%        minimum optode separation in cm.
%   [- invalidate_paired_channels: int, default: 1]
%        When a channel is tagged as bad, also remove the other paired 
%        channels
%   [- nirs_chan_flags: array of bool, default: ones(nb_channels, 1)]
%        Treat only channels where flag is 1. Used to avoid treating
%        auxiliary channels for example.
%  
% Output:
%    - channel_flags: array of int, size: nb_channels
%    - bad_channel_names: cell array of str, size: nb of bad channels

    prev_channel_flags = sData.ChannelFlag;
    channel_flags   = sData.ChannelFlag;
    nirs_flags = strcmpi({channel_def.Channel.Type}, 'NIRS');
    
    signal=sData.F';
    nirs_signal=signal(:,nirs_flags);
    
    nb_chnnels=size(signal,2);
    nb_sample= size(signal,1);

    neg_channels = nirs_flags & any(signal < 0, 1);
    channel_flags(neg_channels) = -1;
    msg(1,:)= {'negative channels', sum(neg_channels)};

    if options.option_sci.Value 
        warning('SCI not available yet');
    end
    
    if options.option_coefficient_variation.Value
        CV_threshold = options.coefficient_variation.Value{1};
        
        CV = std(nirs_signal,1)./mean(nirs_signal,1).*100;
        
        CV_channels = false(1,nb_chnnels);
        CV_channels(nirs_flags) = any(CV>CV_threshold);
        
        channel_flags(CV_channels) = -1;
        msg(end+1,:)= {'high coeffience variation channels', sum(CV_channels)};
    end    
    
    if options.option_remove_saturating.Value
        max_sat_prop = options.option_max_sat_prop.Value{1};
        min_sat_prop = options.option_min_sat_prop.Value{1};
        
        saturating=false(1,nb_chnnels);
        if max_sat_prop < 100
            max_nirs= repmat(max(nirs_signal , [], 1),nb_sample, 1);
            prop_sat_ceil = 100*sum(nirs_signal == max_nirs, 1) / nb_sample;
            saturating(nirs_flags) = saturating(nirs_flags) | prop_sat_ceil >= max_sat_prop;
        end

        if min_sat_prop < 100
            min_nirs= repmat(min(nirs_signal, [], 1), nb_sample , 1);
            prop_sat_floor = 100*sum(nirs_signal == min_nirs, 1) / nb_sample;
            
            saturating(nirs_flags) = saturating(nirs_flags) | prop_sat_floor >= min_sat_prop;
        end
        channel_flags(saturating) = -1;
        msg(end+1,:)= {'saturating channels', sum(saturating)};
    end

    if options.option_separation_filtering.Value
        
        min_separation_m = options.option_separation.Value{1}(1);
        max_separation_m = options.option_separation.Value{1}(2);
        
        separations_m_by_chans = process_nst_separations('Compute', channel_def.Channel(nirs_flags))';
        if all(separations_m_by_chans < 1) % convert to cm if in meter 
            separations_m_by_chans = separations_m_by_chans*100;
        end    
        distances_flag=false(1,nb_chnnels);
        if min_separation_m > 0
            distances_flag(nirs_flags) = distances_flag(nirs_flags) | separations_m_by_chans <= min_separation_m;
        end
        if max_separation_m > 0
             distances_flag(nirs_flags) = distances_flag(nirs_flags) | separations_m_by_chans >= max_separation_m;
        end
        channel_flags(distances_flag) = -1;
        msg(end+1,:)= {'Extrem separation', sum(distances_flag)};
    end 
   
   if ~options.option_keep_unpaired.Value
        channel_flags(nirs_flags) = fix_chan_flags_wrt_pairs(channel_def.Channel(nirs_flags), ...
                                                            channel_def.Nirs.Wavelengths, ...
                                                            channel_flags(nirs_flags) * -1) * -1;
   end
   
   if options.auxilary_signal.Value{1} == 2
       flat_aux =  false(1,nb_chnnels);
       flat_aux(~nirs_flags)= std(signal(:,~nirs_flags),1) < 1; % might need a better threshold 
       
       channel_flags(flat_aux) = -1;
   elseif options.auxilary_signal.Value{1} == 3
       channel_flags(~nirs_flags) = -1;
   end 
   removed =  (prev_channel_flags ~= -1 & channel_flags == -1);
   removed_channel_names = {channel_def.Channel(removed).Name};
   
end

function fixed_chan_flags = fix_chan_flags_wrt_pairs(channel_def, wls, chan_flags)
% Make flags consistent: if flag of a channel is 1, set to 1 all channels
% involved in the same pair

fixed_chan_flags = chan_flags;
nb_wavelengths = length(wls);
nb_channels = length(channel_def);
ichan_to_scan = find(chan_flags==1);
for ii=1:length(ichan_to_scan)
    ichan = ichan_to_scan(ii);
    chan_name = channel_def(ichan).Name;
    pair_prefix = chan_name(1:strfind(chan_name, 'WL'));
    nb_fixed_chans = 0;
    search = {ichan+1:nb_channels ; 1:ichan-1};
    for isearch=1:length(search)
        for i_other_chan=search{isearch}
            if ~isempty(strfind(channel_def(i_other_chan).Name, pair_prefix))
                fixed_chan_flags(i_other_chan) = 1;
                nb_fixed_chans = nb_fixed_chans + 1;
            end
            if nb_fixed_chans == nb_wavelengths-1
                break;
            end
        end
        if nb_fixed_chans == nb_wavelengths-1
            break;
        end
    end
    if nb_fixed_chans ~= nb_wavelengths-1
        throw(MException('NSTError:InconsistentChannel', ...
                         ['Channels paired to ' chan_name ' were not all flagged']));
    end
end
end
