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
% Authors: Edouard Delaire, 2021; Thomas Vincent, 2015-2019


eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    %TOCHECK: how do we limit the input file types (only NIRS data)?
    sProcess.Comment     = 'Detect bad channels';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Pre-process'};
    sProcess.Index       = 1301; %0: not shown, >0: defines place in the list of processes
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
    
    sProcess.options.power_threshold.Comment = 'Power threshold:';
    sProcess.options.power_threshold.Type    = 'value';
    sProcess.options.power_threshold.Value   = {10, '%', 0};
    sProcess.options.power_threshold.Class='sci';

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
    

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Compute bad channel for each Input file
    % Assume that the same mobtage is used for each file
    
    bad_chan_names  = cell(1, length(sInputs));
    bad_chan_msg = cell(1, length(sInputs));
    nirs_flags = cell(1, length(sInputs));

    all_bad_channels = {}; 
    for i_input = 1:length(sInputs)
        % Load channel file
        ChanneMat = in_bst_channel(sInputs(i_input).ChannelFile);
    
        % Load recordings
        if strcmp(sInputs(i_input).FileType, 'data')     % Imported data structure
            sData = in_bst_data(sInputs(i_input).FileName);
        elseif strcmp(sInputs(i_input).FileType, 'raw')  % Continuous data file       
            sData = in_bst(sInputs(i_input).FileName, [], 1, 1, 'no');
        end
        
        nirs_flag= strcmpi({ChanneMat.Channel.Type}, 'NIRS');
        nirs_flags{i_input} = nirs_flag;
        nb_channel=sum(nirs_flag);

        [new_ChannelFlag, bad_chan,msg] = Compute(sData, ChanneMat,sProcess.options);
        bad_chan_names{i_input} = bad_chan;
        bad_chan_msg{i_input} = msg;
        
        all_bad_channels = union(all_bad_channels, bad_chan);
    end
    
    % Set the bad channel in all inputs
    for i_input = 1:length(sInputs)
        tree_set_channelflag({sInputs(i_input).FileName}, 'AddBad', all_bad_channels);
    end
    
    
    % Compute stats  
    % compute the channels that are bad in all runs
    comon_bad_channel = all_bad_channels; 
    for i_input = 1:length(sInputs)
        % Only report NIRS channels
        ChanneMat = in_bst_channel(sInputs(i_input).ChannelFile);
        comon_bad_channel = intersect( comon_bad_channel, intersect( bad_chan_names{i_input}, {ChanneMat.Channel(nirs_flags{i_input}).Name} ));
    end
   
   ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
   all_bad_channels_nirs = intersect( all_bad_channels,  {ChanneMat.Channel(nirs_flags{1}).Name});
   
   fraction= round(100* length(all_bad_channels_nirs)/nb_channel);
   if fraction > 80  % send an error if more than 80% of the channels are marked bad
        bst_report('Error',    sProcess, sInputs, sprintf('%d Channels removed from the files(%d%% of the channels)',length(all_bad_channels_nirs),fraction))
   elseif fraction> 20 % send a warning if more than 20% of the channels are marked bad
        bst_report('Warning',    sProcess, sInputs, sprintf('%d Channels removed from the files(%d%% of the channels)',length(all_bad_channels_nirs),fraction))
   elseif fraction > 0
        bst_report('Info',    sProcess, sInputs, sprintf('%d Channels removed from the files(%d%% of the channels)',length(all_bad_channels_nirs),fraction))
   else 
        bst_report('Info',    sProcess, sInputs, sprintf('No bad channel removed from the files'))
   end
   
   if fraction > 0
         if length(sInputs) > 1 
            bst_report('Info',    sProcess, sInputs, sprintf('Common bad channels: %s', strjoin(comon_bad_channel, ' ,')));
         else
            bst_report('Info',    sProcess, sInputs, sprintf('Bad channels: %s', strjoin(comon_bad_channel, ' ,')));       
         end      
   end
   for i=1: length(sInputs) 
       C = setdiff( intersect(bad_chan_names{i},{ChanneMat.Channel(nirs_flags{1}).Name}),comon_bad_channel); 
       if ~isempty(C) > 0
          bst_report('Info',    sProcess, sInputs(i), sprintf('Bad channels: %s', strjoin(C, ' ,')));
       end
   end     
   
%     critera=zeros( size(msg,1), nb_channel);
%     text=cell(1,size(msg,1));
%     for i = 1:size(msg,1)
%         critera(i,:)= msg{i,2}(nirs_flags);
%         text{i}=msg{i,1};
%     end    
%     confusion=computeConfusion(critera);
%                                             
%     % Warning if too many channels are removed                                         
%     for i=1:size(msg,1)
%         if any(strfind(text{i},'SCI')) || any(strfind(text{i},'power'))
%             continue
%         end
%         fraction= round(100* confusion(i,i)/nb_channel);
%         if fraction > 80  % send an error if more than 80% of the channels are marked bad
%             bst_report('Error',    sProcess, sInputs, sprintf('Number of %s: %d (%d%% of the channels)', text{i},confusion(i,i),fraction))
%         elseif fraction> 20 % send a warning if more than 20% of the channels are marked bad
%             bst_report('Warning',    sProcess, sInputs, sprintf('Number of %s: %d (%d%% of the channels)', text{i},confusion(i,i),fraction))
%         elseif confusion(i,i)> 0
%             bst_report('Info',    sProcess, sInputs, sprintf('Number of %s: %d', text{i},confusion(i,i)))
%         else 
%             bst_report('Info',    sProcess, sInputs, sprintf('No %s', msg{i,1}))
%         end
%     end    
    bst_report('Open', 'current');                                      
    % Add bad channels
    OutputFiles = {sInputs.FileName};
end


%% ===== Compute =====
function [channel_flags, removed_channel_names,criteria] = Compute(sData, channel_def, options)
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
    nb_nirs  = sum(nirs_flags);
    neg_channels = nirs_flags & any(signal < 0, 1);
    channel_flags(neg_channels) = -1;
    criteria(1,:)= {'negative channels', neg_channels,{} };
    

    if options.option_sci.Value 
        
            fs = 1 / diff(sData.Time(1:2));
            [SCI, power] = process_nst_sci('compute',nirs_signal',channel_def.Channel(nirs_flags), fs);
            
            SCI = SCI *100;
            power = power*100;
            
            SCI_threshold = options.sci_threshold.Value{1};
            SCI_channels = false(1,nb_chnnels);
            SCI_channels(nirs_flags) = SCI < SCI_threshold ;
            
            power_threshold =  options.power_threshold.Value{1};
            power_channels = false(1,nb_chnnels);
            power_channels(nirs_flags) = power < power_threshold ;
            
            %create figure
            fig1= figure('visible','off');
            h1_axes=subplot(1,1,1);
            h=histogram(SCI');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            line(h1_axes,[SCI_threshold, SCI_threshold], ylim(gca), 'LineWidth', 2, 'Color', 'b');
        
            leg=legend(h1_axes,'SCI threshold');
            title(sprintf('SCI rejection: %d channels',sum(SCI_channels))) 
            criteria(end+1,:)= {'Low SCI', SCI_channels,{h1_axes,leg}};

            fig1= figure('visible','off');
            h1_axes=subplot(1,1,1);
            h=histogram(power');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
            line(h1_axes,[power_threshold, power_threshold], ylim(gca), 'LineWidth', 2, 'Color', 'b');
        
            leg=legend(h1_axes,'SCI threshold');
            title(sprintf('power rejection: %d channels',sum(power_channels))) 
            
            criteria(end+1,:)= {'Low power', power_channels,{h1_axes,leg}};
            criteria(end+1,:)= {'Cardiac', SCI_channels & power_channels,{}};

            channel_flags(SCI_channels & power_channels) = -1;
    end
    
    if options.option_coefficient_variation.Value
        CV_threshold = options.coefficient_variation.Value{1};                 
        low_cutoff  = 1/200;
        high_cutoff = 0.1;
        order       = 2;
        fs = 1 / diff(sData.Time(1:2));

        [nirs_filtered, transient] = process_nst_iir_filter('Compute',nirs_signal(100:end,:), fs, 'bandpass', low_cutoff, ...
                                  high_cutoff, order, 1);
        
        %CV2 = std(nirs_signal,1)./mean(nirs_signal,1).*100;
        CV = std(nirs_filtered,1,'omitnan')./mean(nirs_filtered,1).*100;
        
        % To see the impact of adding the filter: (large drop of CV)
        % figure; boxplot([CV,CV2],[ ones(1,size(nirs_signal,2)),2*ones(1,size(nirs_signal,2))])
        
        CV_channels = false(1,nb_chnnels);
        CV_channels(nirs_flags) = CV>CV_threshold ;
        
        channel_flags(CV_channels) = -1;
        
        %create figure
%         fig1= figure('visible','off');
%         h1_axes=subplot(1,1,1);
%         h=histfit(CV',12,'normal');
%         h(1).Annotation.LegendInformation.IconDisplayStyle = 'off';
%         h(2).Annotation.LegendInformation.IconDisplayStyle = 'off';
% %         
%         line(h1_axes,[CV_threshold, CV_threshold], ylim(gca), 'LineWidth', 2, 'Color', 'b');
%         line(h1_axes,[icdf(CV_fit,0.99), icdf(CV_fit,0.99)], ylim(gca), 'LineWidth', 2, 'Color', 'g');
%         
%         leg=legend(h1_axes,{sprintf('Used threshold: %.1f (Prob of rejection %.3f%%)',CV_threshold,100*cdf(CV_fit,CV_threshold,'upper')), ...
%                               sprintf('Suggested threshold: %.3f',icdf(CV_fit,0.99))});
%                            
%         title(sprintf('CV rejection: %d channels',sum(CV_channels)))        
        criteria(end+1,:)= {'high coeffience variation channels', CV_channels,{}};
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
        criteria(end+1,:)= {'saturating channels', saturating,{} };
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
        
        %create figure

        fig1= figure('visible','off');
        h1_axes=subplot(1,1,1);
        h=histogram(separations_m_by_chans');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        line(h1_axes,[min_separation_m, min_separation_m], ylim(gca), 'LineWidth', 2, 'Color', 'b');
        line(h1_axes,[max_separation_m, max_separation_m], ylim(gca), 'LineWidth', 2, 'Color', 'b');
        
        leg=legend(h1_axes,'Separation threshold');
        title(sprintf('Distance rejection: %d channels',sum(distances_flag))) 
        
        
        criteria(end+1,:)= {'Extrem separation', distances_flag,{h1_axes,leg} };
    end 
    
    % plot SCI and power according to distance
    if 0
    figure; 
    h1_axes= subplot(2,1,1); hold on;
    plot(separations_m_by_chans,SCI,'+')
    line(h1_axes,[min_separation_m, min_separation_m], ylim(gca), 'LineWidth', 1, 'Color', 'b');
    line(h1_axes,[max_separation_m, max_separation_m], ylim(gca), 'LineWidth', 1, 'Color', 'b');
    line(h1_axes,xlim(gca), [SCI_threshold, SCI_threshold], 'LineWidth', 2, 'Color', 'r');
    plot(separations_m_by_chans(SCI_channels & power_channels),SCI(SCI_channels & power_channels),'o')
    xlabel('Source-Detecteur separation (cm)');ylabel('SCI(%)');
    title('Rejection based on SCI and power')
    h2_axes=subplot(2,1,2); hold on;
    plot(separations_m_by_chans,power,'+')
    line(h2_axes,[min_separation_m, min_separation_m], ylim(gca), 'LineWidth', 1, 'Color', 'b');
    line(h2_axes,[max_separation_m, max_separation_m], ylim(gca), 'LineWidth', 1, 'Color', 'b');
    line(h2_axes,xlim(gca), [power_threshold, power_threshold], 'LineWidth', 2, 'Color', 'r');   
    plot(separations_m_by_chans(SCI_channels & power_channels),power(SCI_channels & power_channels),'o')
    xlabel('Source-Detecteur separation (cm)');ylabel('power(%)');
    
    figure; % relation CV/Distance
    h1_axes= subplot(1,1,1); hold on;
    plot(separations_m_by_chans,CV,'+')
    line(h1_axes,[min_separation_m, min_separation_m], ylim(gca), 'LineWidth', 1, 'Color', 'b');
    line(h1_axes,[max_separation_m, max_separation_m], ylim(gca), 'LineWidth', 1, 'Color', 'b');
    line(h1_axes,xlim(gca), [CV_threshold, CV_threshold], 'LineWidth', 2, 'Color', 'r');
    plot(separations_m_by_chans(CV > CV_threshold),CV(CV > CV_threshold),'o')
    xlabel('Source-Detecteur separation (cm)');ylabel('CV(%)');
    end

   if ~options.option_keep_unpaired.Value
        channel_flags(nirs_flags) = fix_chan_flags_wrt_pairs(channel_def.Channel(nirs_flags), ...
                                                            channel_def.Nirs.Wavelengths, ...
                                                            channel_flags(nirs_flags) * -1) * -1;
   end
   
   % Remove auxilary signal
   if options.auxilary_signal.Value{1} == 2 % remove flat signal
       flat_aux =  false(1,nb_chnnels);
       flat_aux(~nirs_flags)= std(signal(:,~nirs_flags),1) < 1; % might need a better threshold 
       
       channel_flags(flat_aux) = -1;
   elseif options.auxilary_signal.Value{1} == 3 % remove all signal
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

function confusion=computeConfusion(criteria)
    % Compute the 'confusion' matrix. confusion(i,j) contains the number of
    % channel marked as bad with criteria i, and j (exemple negative
    % channel and channel with high variance 
    % Can be plotted usin confusionchart(m,classLabels) from machine
    % learning toolbox
    
    n_criteria=size(criteria,1);
    
    confusion=zeros(n_criteria,n_criteria);
    for i =1:n_criteria
        for j=1:n_criteria
            confusion(i,j) = sum( criteria(i,:) & criteria(j,:) ); 
        end
    end    

end

function plotConfMat(varargin)
%PLOTCONFMAT plots the confusion matrix with colorscale, absolute numbers
%   and precision normalized percentages
%
%   usage: 
%   PLOTCONFMAT(confmat) plots the confmat with integers 1 to n as class labels
%   PLOTCONFMAT(confmat, labels) plots the confmat with the specified labels
%
%   Vahe Tshitoyan
%   20/08/2017
%
%   Arguments
%   confmat:            a square confusion matrix
%   labels (optional):  vector of class labels
% number of arguments
switch (nargin)
    case 0
       confmat = 1;
       labels = {'1'};
    case 1
       confmat = varargin{1};
       labels = 1:size(confmat, 1);
    otherwise
       confmat = varargin{1};
       labels = varargin{2};
end

confmat(isnan(confmat))=0; % in case there are NaN elements
numlabels = size(confmat, 1); % number of labels

selection = any( confmat > 0, 1) & ~contains(labels,'SCI') & ~contains(labels,'power');


% calculate the percentage accuracies
tot_rejection = sum( diag( confmat(selection,selection))); 
confpercent = 100*confmat ./tot_rejection;

selection = any( confmat > 0, 1) & ~contains(labels,'Cardiac');
if all(confmat(contains(labels,'power'),:) == 0)  &&  any(confmat(contains(labels,'SCI'),:)) > 0
    selection( find(contains(labels,'power'))) = 1;
elseif all(confmat(contains(labels,'SCI'),:) == 0)  &&  any(confmat(contains(labels,'power'),:)) > 0
    selection( find(contains(labels,'SCI'))) = 1;
end    
confpercent=confpercent(selection,selection);
confmat = confmat(selection,selection);
labels = labels(selection);
numlabels= length(labels);

% plotting the colors
imagesc(confpercent);
title(sprintf('Total rejected channel: %d', tot_rejection));
ylabel('Rejection criteria'); xlabel('Rejection criteria');
% set the colormap
colormap(flipud(gray));
% Create strings from the matrix values and remove spaces 
textStrings = num2str([confpercent(:), confmat(:)], '%.1f%%\n%d\n');

textStrings = strtrim(cellstr(textStrings));
% Create x and y coordinates for the strings and plot them
[x,y] = meshgrid(1:numlabels);
hStrings = text(x(:),y(:),textStrings(:), ...
    'HorizontalAlignment','center');
% Get the middle value of the color range
midValue = mean(get(gca,'CLim'));
% Choose white or black for the text color of the strings so
% they can be easily seen over the background color
textColors = repmat(confpercent(:) > midValue,1,3);
set(hStrings,{'Color'},num2cell(textColors,2));
% Setting the axis labels
set(gca,'XTick',1:numlabels,...
    'XTickLabel',labels,...
    'YTick',1:numlabels,...
    'YTickLabel',labels,...
    'TickLength',[0 0]);
end

