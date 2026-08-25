function varargout = process_nst_bootstrap_MNE( varargin )
% process_nst_bootstrap_MNE: use bootstrap estimator of the mean followed by MNE.
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2018 University of Southern California & McGill University
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
% Authors: Edouard Delaire 2026, Jean-Eudes Bornert 2026
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Boostrap MNE';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', '3D reconstruction'};
    sProcess.Index       = 1505;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    
    sProcess.options.label0.Comment = '<b>Bootstraps options:</b>';
    sProcess.options.label0.Type = 'label';

    
    sProcess.options.timewindow_snr.Comment = 'Time window for SNR:';
    sProcess.options.timewindow_snr.Type    = 'timewindow';
    sProcess.options.timewindow_snr.Value   =  [];
    sProcess.options.timewindow_snr.Group   = 'MNE';

    sProcess.options.replacement.Comment = 'Do bootstrap with replacement ?';
    sProcess.options.replacement.Type    = 'checkbox';
    sProcess.options.replacement.Value   = 1;

    sProcess.options.combination.Comment = 'Combination value #:';
    sProcess.options.combination.Type    = 'value';
    sProcess.options.combination.Value   = {16, '', 0};

    sProcess.options.n_average.Comment = 'Number of average:';
    sProcess.options.n_average.Type    = 'value';
    sProcess.options.n_average.Value   = {50, '', 0};

    sProcess.options.label1.Comment = '<b>MNE options:</b>';
    sProcess.options.label1.Type = 'label';

    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {5, 'cm',2};
    sProcess.options.thresh_dis2cortex.Group   = 'MNE';

    sProcess.options.depth_weightingMNE.Comment = str_pad('Depth weighting factor for MNE',35);
    sProcess.options.depth_weightingMNE.Type    = 'value';
    sProcess.options.depth_weightingMNE.Value   = {0.5, '', 1};
    sProcess.options.depth_weightingMNE.Group   = 'MNE';

    sProcess.options.TimeSegment.Comment = str_pad('Reconstruction Time window:',35);
    sProcess.options.TimeSegment.Type    = 'timewindow';
    sProcess.options.TimeSegment.Value   = [];
    sProcess.options.TimeSegment.Group   = 'MNE';

    sProcess.options.NoiseCov_recompute.Comment = 'Compute noise covariance of the baseline MNE';
    sProcess.options.NoiseCov_recompute.Type    = 'checkbox';
    sProcess.options.NoiseCov_recompute.Controller = 'noise_cov';
    sProcess.options.NoiseCov_recompute.Value   = 1;
    sProcess.options.NoiseCov_recompute.Group   = 'MNE';

    sProcess.options.TimeSegmentNoise.Comment = str_pad('Baseline Time window:',35);
    sProcess.options.TimeSegmentNoise.Type    = 'baseline';
    sProcess.options.TimeSegmentNoise.Value   = [];
    sProcess.options.TimeSegmentNoise.Class = 'noise_cov';
    sProcess.options.TimeSegmentNoise.Group   = 'MNE';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end

function s = str_pad(s,padsize)
    if (length(s) < padsize)
        s = [repmat('&nbsp;', 1, padsize - length(s)), s];
    end
    s = ['<FONT FACE="monospace">' s '</FONT>'];
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    
    % Load data
    sStudy = bst_get('Study', sInputs.iStudy);
    if isempty(sStudy.iHeadModel)
        bst_error('No head model found. Consider running "NIRS -> Compute head model"');
        return;
    end

    % Load Head model
    HeadModelFileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    sHead = in_bst_headmodel(HeadModelFileName, 1);

    % Load Data
    [sChannel, dataMat, Time] = LoadData(sInputs);

    % Prepare options
    nTrials = length(sInputs);
    options = getOptions(sProcess, Time, nTrials);

    % Generate all permutations; and compute SNR for each permutation
    [avg_list, SNR] = generate_permutations(sChannel, dataMat,  options);
    
    % Select the permutation around the median SNR for all wavelenght (e.g. 101 resamples)
    [avg_list, SNR_selected, quantiles] = generate_avg_list(avg_list, SNR, options.nAverage);
    plot_trials(avg_list, SNR, SNR_selected, quantiles);
    
    % Genrate median average
    sProcess.options.avg_func.Value = 7; % Arithmetic average + Standard error
    sFiles = process_average('AverageFiles', sProcess, sInputs(avg_list(1, :)), 1);

    % sFiles = bst_process('CallProcess', 'process_average', sInputs(avg_list(1, :)), [], ...
    %     'avgtype',       5, ...  % By trial group (folder average)
    %     'avg_func',      7, ...  
    %     'weighted',      0, ...
    %     'keepevents',    1);
    
    % Process: Add tag:  
    AvgFile = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
        'tag',           sprintf( '| median SNR = %.2f', median(SNR_selected(1))) , ...
        'output',        1);  % Add to file name
    
    
    sDataIn = in_bst_data(AvgFile.FileName);
    ChannelMat = in_bst_channel(AvgFile.ChannelFile);
    OPTIONS = process_nst_wmne('getOptions', sProcess, HeadModelFileName, AvgFile.FileName);
    
        
    bst_progress('start', 'Bootstraping', 'Reconstruction by MNE', 1, size(avg_list, 1)); 
    Results = zeros(size(avg_list, 1) ,size(sHead.Gain, 2), 2, length(sDataIn.Time));
    for iAvg = 1:size(avg_list, 1)
        
        avg_trial = mean(dataMat(avg_list(iAvg,:), : , :), 1);
        sDataIn.F = squeeze(avg_trial);
    
        [sResults] = process_nst_wmne('Compute', OPTIONS, ChannelMat, sDataIn);
        sResults = process_nst_wmne('filterResults', sResults, [0, 1, 1, 0]);
    
        Results(iAvg,:,1,:) = bst_multiply_cellmat(sResults(1).ImageGridAmp);
        Results(iAvg,:,2,:) = bst_multiply_cellmat(sResults(2).ImageGridAmp);
        
        bst_progress('inc', 1);
    end 
    
    ResultsAvg = squeeze(mean(Results, 1));
    ResultsSD  = squeeze(std(Results, [], 1));
    
    % Save results
    bst_progress('text', 'Saving Results...');
    for iMap = 1:length(sResults)
        ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' nst_protect_fn_str(sResults(iMap).Comment)]);
    
        ResultsMat          = sResults(iMap);
        ResultsMat.DataFile = AvgFile.FileName;
        ResultsMat.Options  = OPTIONS;
        
        ResultsMat.ImageGridAmp = squeeze(ResultsAvg(:,iMap,:));
        ResultsMat.Std = squeeze(ResultsSD(:,iMap,:));
    
        ResultsMat.nAvg = options.nAverage;
        ResultsMat.Leff = options.nAverage;
    
        bst_save(ResultFile, ResultsMat, 'v6');
        db_add_data(sInputs(1).iStudy, ResultFile, ResultsMat);
    
        OutputFiles{end+1} = ResultFile;
    end
end

function plot_trials(avg_list, SNR, SNR_selected, quantiles)

    % TODO: replace hist with histogram

    figure
    subplot(2,2,1);
    histogram(avg_list,'BinMethod','integers')
    % line(xlim(gca), length(avg_all_list)*(sProcess.options.combination.Value{1}/length(sInputs))*[1,1],'Color','red','LineStyle','--')
    xlabel('Trials');
    title('Occurance of each trial in the final result')
    
    subplot(2, 2, 2);
    hist(SNR(1, :), size(SNR, 2)./10);
    line(median(SNR_selected(1,:))*[1,1],ylim(gca),'Color','red','LineStyle','--')
    line(min(SNR_selected(1,:))*[1,1],ylim(gca),'Color','blue','LineStyle','--')
    line(max(SNR_selected(1,:))*[1,1],ylim(gca),'Color','blue','LineStyle','--')
    xlim([0 43])
    title('SNR for \lambda = 690');
    
    subplot(2,2,3);
    hist(SNR(2, :), size(SNR,2)./10);
    line(median(SNR_selected(2,:))*[1,1],ylim(gca),'Color','red','LineStyle','--')
    line(min((SNR_selected(2,:)))*[1,1],ylim(gca),'Color','blue','LineStyle','--')
    line(max((SNR_selected(2,:)))*[1,1],ylim(gca),'Color','blue','LineStyle','--')
    xlim([0 40])
    title('SNR for \lambda = 830');
    
    subplot(2,2,4);
    hist(mean(SNR, 1), size(SNR,2)./10);
    line(median(mean(SNR_selected))*[1,1],ylim(gca),'Color','red','LineStyle','--')
    line(min(mean(SNR_selected))*[1,1],ylim(gca),'Color','blue','LineStyle','--')
    line(max(mean(SNR_selected))*[1,1],ylim(gca),'Color','blue','LineStyle','--')
    xlim([0 40])
    title('Average SNR for \lambda = {690, 830}');
    
    figure;
    scatterhist(SNR(1, :), SNR(2, :),'Location','NorthEast', 'Direction','out'); hold on;
    scatter(SNR_selected(1, :), SNR_selected(2, :),'filled','r')
    yline(quantiles(2, :) ); xline(quantiles(1, :))
    sgtitle('SNR for \lambda = {690, 830}')
end

function [avg_list, SNR] = generate_permutations(ChannelMat, DataMat_dOD, options)

    nTrials = size(DataMat_dOD, 1);

    if options.isExact
        avg_list    = nchoosek(1:nTrials, options.nCombination); 
        n_perm      = size(avg_list, 1);
    else
        n_perm      = 5000;
        avg_list    = zeros(n_perm, options.nCombination);
    end
    
    bst_progress('start', 'Bootstraping', 'Computing all combinations', 0, length(avg_list)); 
    
    SNR = zeros(2, n_perm);
    
    for iAvg = 1:n_perm
    
        if ~options.isExact
            avg_list(iAvg,:) =  randsample(nTrials, options.nCombination, options.isReplacement); 
        end
    
        avg_trial = squeeze(mean(DataMat_dOD(avg_list(iAvg,:), : , :), 1));
    
        trials_690 = avg_trial(1:2:end,:);
        trials_830 = avg_trial(2:2:end,:);
    
        trials_std_690 = std(trials_690(:,ibaseline),[],2);
        trials_std_830 = std(trials_830(:,ibaseline),[],2);
    
        SNR(1, iAvg) = max(max(abs(trials_690(:, iSNR)))./mean(trials_std_690));
        SNR(2, iAvg) = max(max(abs(trials_830(:, iSNR)))./mean(trials_std_830));

        bst_progress('inc', 1);
    end

    bst_progress('stop');
end

function [avg_list, SNR, quantiles] = generate_avg_list(avg_list, SNR, nAverage)

    nAverage = min(nAverage, floor(length(avg_list)/2-1));
    
    quantiles_list = [0.5 0.6 0.7 0.8 0.9];
    quantiles = zeros(size(SNR, 1), length(quantiles_list));

    for iWavelenght = 1:size(SNR, 1)
        quantiles(iWavelenght, :) = quantile(SNR(iWavelenght,:)', quantiles_list);
    end

    %this is targeting median SNR for both wavelength
    SNR_all = sqrt(sum( (SNR - repmat(quantiles(:, 1), 1, size(SNR, 2))).^2 , 1));

    
    [SNR_sorted, SNR_all_list] = sort(SNR_all);
    [~, idx_list] =  min( abs(SNR_sorted - min(SNR_sorted)) );
    
    idx_selected = SNR_all_list(idx_list:min(length(SNR_all), idx_list+2*nAverage));

    avg_list = avg_list(idx_selected, :);
    SNR  = SNR(:, idx_selected);
    
end

function [sChannel, dataMat, Time] = LoadData(sInputs)

    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    
    nTrials = length(sInputs);
    sData   = cell(1, nTrials);
    for iFile = 1:nTrials
        sData{iFile} = in_bst_data(sInputs(iFile).FileName, 'F', 'ChannelFlag', 'History', 'Time');
    end
    
    iChannel    = good_channel(ChannelMat.Channel, sData{1,1}.ChannelFlag, 'NIRS'); 
    Time        = sData{1,1}.Time;
    
    dataMat = zeros(nTrials, length(iChannel), length(Time));
    sChannel = ChannelMat.Channel(iChannel);
    
    for iFile = 1:nTrials
        
        if ~all(sData{iFile}.ChannelFlag == sData{1}.ChannelFlag)
            bst_error('All trials must have the same bad channels');
        end
        
        dataMat(iFile, :, :) = sData{iFile}.F(iChannel, :);
    end

end

function options = getOptions(sProcess, Time, nTrials)

    timewindow_baseline = sProcess.options.TimeSegmentNoise.Value{1};
    timewindow_snr      = sProcess.options.timewindow_snr.Value{1};
    isReplacement = sProcess.options.replacement.Value;
    nCombination  = sProcess.options.combination.Value{1};
    nAverage      = sProcess.options.n_average.Value{1};
    
    
    ibaseline = panel_time('GetTimeIndices', Time, timewindow_baseline);
    iSNR      = panel_time('GetTimeIndices', Time, timewindow_snr);
    isExact   = nTrials <= 20  &&  ~isReplacement;
    options = struct('nCombination', nCombination, 'nAverage', nAverage, 'isReplacement', isReplacement, 'isExact', isExact, 'ibaseline', ibaseline, 'iSNR', iSNR);
end