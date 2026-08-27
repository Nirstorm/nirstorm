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
    sProcess.options.timewindow_snr.Type    = 'poststim';
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


    sProcess.options.target_quantile.Comment = {'median', '70 percentile', '80 percientile', '90 percentile', 'targeted SNR quantile:'; '0.5', '0.7', '0.8', '0.9',  ''} ;
    sProcess.options.target_quantile.Type    = 'radio_linelabel';
    sProcess.options.target_quantile.Value   = '0.5';

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

    sProcess.options.output_all.Comment = 'Save all bootstrap samples (warning: requires a lot of memory)';
    sProcess.options.output_all.Type    = 'checkbox';
    sProcess.options.output_all.Value   = 1;
    sProcess.options.output_all.Group   = 'output';

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
    
    
    % Choose appropriate inverse function
    if isfield(sProcess.options, 'mem') && isfield(sProcess.options.mem, 'Value') && ~isempty(sProcess.options.mem.Value)
        inverse_function = @process_nst_cmem;
        function_name    = 'cMEM';
    else
        inverse_function = @process_nst_wmne;
        function_name    = 'MNE';
    end


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
    [ChannelMat, sChannel, dataMat, Time] = LoadData(sInputs);

    % Prepare options
    nTrials = length(sInputs);
    options = getOptions(sProcess, Time, nTrials);

    % Generate all permutations; and compute SNR for each permutation
    [avg_list, SNR] = generate_permutations(sChannel, dataMat,  options);
    
    % Select the permutation around the median SNR for all wavelenght (e.g. 101 resamples)
    [avg_list, SNR_selected, quantiles] = generate_avg_list(avg_list, SNR, options.nAverage, options.target_snr);
    plot_trials(avg_list, SNR, SNR_selected, quantiles);
    
    % Genrate median average
    AvgFile = computeAverage(sInputs(avg_list(1, :)), sprintf( 'median SNR = %.2f / %.2f', median(SNR_selected(1, :)) ,  median(SNR_selected(2, :))));

    sDataIn = in_bst_data(AvgFile.FileName);
    OPTIONS = inverse_function('getOptions', sProcess, HeadModelFileName, AvgFile.FileName);

    bst_progress('start', 'Bootstraping', sprintf('Reconstruction by %s', function_name), 1, size(avg_list, 1)); 

    meanvar_estimator(1) = nst_math_WelfordVariance(size(avg_list, 1), [size(sHead.Gain, 2), length(sDataIn.Time)]);
    meanvar_estimator(2) = nst_math_WelfordVariance(size(avg_list, 1), [size(sHead.Gain, 2), length(sDataIn.Time)]);

    for iAvg = 1:size(avg_list, 1)
        % Put the data in place
        sDataIn.F   = squeeze(mean(dataMat(avg_list(iAvg,:), : , :), 1));
        sDataIn.Std = []; 

        % Compute inverse solution (MNE or cMEM)
        sResults = inverse_function('Compute', OPTIONS, ChannelMat, sDataIn);
        % Select HbO and HbR
        sResults = inverse_function('filterResults', sResults, [0, 1, 1, 0]); 
        
        % Store resuls
        meanvar_estimator(1).update(bst_multiply_cellmat(sResults(1).ImageGridAmp))
        meanvar_estimator(2).update(bst_multiply_cellmat(sResults(2).ImageGridAmp))
        
        if options.output_all
            tmp = computeAverage(sInputs(avg_list(iAvg, :)), sprintf('bootstrap (#%d)', iAvg));
            for iMap = 1:length(sResults)

                ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' nst_protect_fn_str(sResults(iMap).Comment)]);
            
                ResultsMat          = sResults(iMap);
                ResultsMat.DataFile = file_short(tmp.FileName);
                ResultsMat.Options  = OPTIONS;
                
                ResultsMat.nAvg = options.nCombination;
                ResultsMat.Leff = options.nCombination;
            
                bst_save(ResultFile, ResultsMat, 'v6');
                db_add_data(sInputs(1).iStudy, ResultFile, ResultsMat);
            
                OutputFiles{end+1} = ResultFile;
            end
        end
    
        bst_progress('inc', 1);
    end 
    
    % Save results
    bst_progress('text', 'Saving Results...');
    for iMap = 1:length(sResults)
        ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' nst_protect_fn_str(sResults(iMap).Comment)]);
    
        ResultsMat          = sResults(iMap);
        ResultsMat.DataFile = file_short(AvgFile.FileName);
        ResultsMat.Options  = OPTIONS;
        
        ResultsMat.ImageGridAmp = meanvar_estimator(iMap).Mean;
        ResultsMat.Std = meanvar_estimator(iMap).StdDev;
    
        ResultsMat.nAvg = meanvar_estimator(iMap).N;
        ResultsMat.Leff = meanvar_estimator(iMap).N;
    
        bst_save(ResultFile, ResultsMat, 'v6');
        db_add_data(sInputs(1).iStudy, ResultFile, ResultsMat);
    
        OutputFiles{end+1} = ResultFile;
    end
end

function plot_trials(avg_list, SNR, SNR_selected, quantiles)

    p = figure;
    u1 = uipanel(p, 'position', [0, 0, 0.5, 1]);
    ax1 = axes('Parent', u1, 'Position', [0.15, 0.15, 0.75, 0.75]);

    histogram(ax1, avg_list,'BinMethod','integers');
    xlabel('Trials');
    title('Occurance of each trial in the final result')
    

    u2 = uipanel(p, 'position', [0.5, 0, 0.5, 1]);  
    h = scatterhist(SNR(1, :), SNR(2, :),'Location','NorthEast', 'Direction','out', 'Parent', u2);

    hold(h(1), 'on'); 
    scatter(h(1), SNR_selected(1, :), SNR_selected(2, :), 'filled', 'r');
    yline(h(1), quantiles(2, :) ); xline(h(1), quantiles(1, :))
    sgtitle('SNR for \lambda')

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
    
    ibaseline = options.ibaseline;
    iSNR      = options.iSNR;

    bst_progress('start', 'Bootstraping', 'Computing all combinations', 0, length(avg_list)); 

    groups = unique({ChannelMat.Group});
    SNR = zeros(length(groups), n_perm);
    
    for iAvg = 1:n_perm
        if ~options.isExact
            avg_list(iAvg,:) =  randsample(nTrials, options.nCombination, options.isReplacement); 
        end
    
        avg_trial = squeeze(mean(DataMat_dOD(avg_list(iAvg,:), : , :), 1));
        for iGroup = 1:length(groups)
            iChannels = strcmp({ChannelMat.Group}, groups{iGroup});
            
            avg_channel     = avg_trial(iChannels, :);
            snr_max         = max(abs(avg_channel(:, iSNR)), [], 2); 
            baseline_std    = std(avg_channel(:,ibaseline), [], 2);

            SNR(iGroup, iAvg) = max(snr_max) ./ mean(baseline_std);
        end
        bst_progress('inc', 1);
    end

    bst_progress('stop');
end

function [avg_list, SNR, quantiles] = generate_avg_list(avg_list, SNR, nAverage, targeted_snr)

    nAverage = min(nAverage, floor(length(avg_list)/2-1));
    
    quantiles_list = [0.5 0.6 0.7 0.8 0.9];
    quantiles = zeros(size(SNR, 1), length(quantiles_list));

    for iWavelenght = 1:size(SNR, 1)
        quantiles(iWavelenght, :) = quantile(SNR(iWavelenght,:)', quantiles_list);
    end

    % Select sample according to the targeted SNR 
    [~, iQuantile] = min( abs(quantiles_list - targeted_snr));
    SNR_all = sqrt(sum( (SNR - repmat(quantiles(:, iQuantile), 1, size(SNR, 2))).^2 , 1));

    
    [SNR_sorted, SNR_all_list] = sort(SNR_all);
    [~, idx_start]  = min(abs(SNR_sorted - min(SNR_sorted)));
    idx_end         = min(length(SNR_all), idx_start + nAverage);
    idx_selected    = SNR_all_list(idx_start:idx_end);

    avg_list = avg_list(idx_selected, :);
    SNR      = SNR(:, idx_selected);
    
end

function [ChannelMat, sChannel, dataMat, Time] = LoadData(sInputs)

    ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
    
    nTrials = length(sInputs);
    sData   = cell(1, nTrials);
    for iFile = 1:nTrials
        sData{iFile} = in_bst_data(sInputs(iFile).FileName, 'F', 'ChannelFlag', 'History', 'Time');
    end
    
    iChannel    = good_channel(ChannelMat.Channel, sData{1,1}.ChannelFlag, 'NIRS'); 
    nChannel    = length(sData{1,1}.ChannelFlag);
    Time        = sData{1,1}.Time;
    

    dataMat = zeros(nTrials, nChannel, length(Time));
    sChannel = ChannelMat.Channel;

    for iFile = 1:nTrials
        
        if ~all(sData{iFile}.ChannelFlag == sData{1}.ChannelFlag)
            bst_error('All trials must have the same bad channels');
        end
        
        data_trial = nan(nChannel, length(Time));
        data_trial(iChannel, :) = sData{iFile}.F(iChannel, :);

        dataMat(iFile, :, :) = data_trial;
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
    
    output_all = sProcess.options.output_all.Value;

    target_snr = str2double(sProcess.options.target_quantile.Value);
    options = struct('nCombination', nCombination, 'nAverage', nAverage, ...
                      'isReplacement', isReplacement, 'isExact', isExact, ...
                      'ibaseline', ibaseline, 'iSNR', iSNR, ...
                      'target_snr', target_snr, ...
                      'output_all', output_all);
end

function AvgFile = computeAverage(sInputs, Tag)
% computeAverage; Compute the average of the sInputs file, and add Tag to
% the file comment. 

    % Process: Compute average
    sProcess =  process_average('GetDescription');
    sProcess.options.avg_func.Value = 7; % Arithmetic average + Standard error
    sProcess.options.avgtype.Value  = 1;
    sFiles = process_average('AverageFiles', sProcess, sInputs, 1);
    
    % Prepare input of add tag.
    sInputs = sInputs(1);
    sInputs.FileName = sFiles;
    
    if isempty(Tag)
        % If there is not tag to add. 
        AvgFile = sInputs;
        AvgFile.FileName = file_short(sFiles);

        return;
    end

    % Process: Add tag:  
    sProcess = process_add_tag('GetDescription');
    sProcess.options.tag.Value = Tag; 
    sProcess.options.output.Value = 'name';
    tmp = process_add_tag('Run', sProcess, sInputs);

    
    % Prepare output
    AvgFile = sInputs;
    AvgFile.FileName = file_short(tmp{1});
end