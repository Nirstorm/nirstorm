function varargout = process_nst_bootstrap_MEM( varargin )
% process_nst_bootstrap_MNE:  Compute source-localization using Bootstrap
% and MEM

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
% Authors: Edouard Delaire, 2022

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Boostrap MEM';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', '3D reconstruction'};
    sProcess.Index       = 1506;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'ressults'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;


    sProcess.options.label0.Comment = '<b>Bootstraps options:</b>';
    sProcess.options.label0.Type = 'label';

    sProcess.options.TimeSegmentNoise.Comment = 'Baseline Time window:';
    sProcess.options.TimeSegmentNoise.Type    = 'baseline';
    sProcess.options.TimeSegmentNoise.Value   = [];

    sProcess.options.timewindow_snr.Comment = 'Time window for SNR:';
    sProcess.options.timewindow_snr.Type    = 'timewindow';
    sProcess.options.timewindow_snr.Value   =  [];

    sProcess.options.replacement.Comment = 'Do bootstrap with replacement ?';
    sProcess.options.replacement.Type    = 'checkbox';
    sProcess.options.replacement.Value   = 1;

    sProcess.options.combination.Comment = 'Combination value #:';
    sProcess.options.combination.Type    = 'value';
    sProcess.options.combination.Value   = {16, '', 0};

    sProcess.options.n_average.Comment = 'Number of average:';
    sProcess.options.n_average.Type    = 'value';
    sProcess.options.n_average.Value   = {50, '', 0};

    % Time window for baseline
    sProcess.options.combination.Comment = 'Combination value #:';
    sProcess.options.combination.Type    = 'value';
    sProcess.options.combination.Value   = {16, '', 0};
    
    % Definition of the options for source reconstruction
    sProcess.options.label1.Comment = '<b>MEM options:</b>';
    sProcess.options.label1.Type = 'label';

    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {5, 'cm',2};
    sProcess.options.thresh_dis2cortex.Group   = 'MNE';


    sProcess.options.mem.Comment = {'panel_brainentropy', 'Source estimation options: '};
    sProcess.options.mem.Type    = 'editpref';
    sProcess.options.mem.Value   = be_main();

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};

    % Install/load brainentropy plugin
    [isInstalled, errMessage] = bst_plugin('Install', 'brainentropy', 1);
    if ~isInstalled
        return;
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
    [ChannelMat, sChannel, dataMat, Time] = process_nst_bootstrap_MNE('LoadData', sInputs);


    % Prepare options
    nTrials = length(sInputs);
    options = process_nst_bootstrap_MNE('getOptions', sProcess, Time, nTrials);


    % Generate all permutations; and compute SNR for each permutation
    [avg_list, SNR] = process_nst_bootstrap_MNE('generate_permutations', sChannel, dataMat,  options);
    
    % Select the permutation around the median SNR for all wavelenght (e.g. 101 resamples)
    [avg_list, SNR_selected, quantiles] = process_nst_bootstrap_MNE('generate_avg_list', avg_list, SNR, options.nAverage);
    process_nst_bootstrap_MNE('plot_trials', avg_list, SNR, SNR_selected, quantiles);

    % Genrate median average
    AvgFile = process_nst_bootstrap_MNE('computeAverage', sInputs(avg_list(1, :)), sprintf( 'median SNR = %.2f', median(SNR_selected(1))));

    sDataIn = in_bst_data(AvgFile.FileName);
    OPTIONS = process_nst_cmem('getOptions', sProcess, HeadModelFileName, AvgFile.FileName);
    
    bst_progress('start', 'Bootstraping', 'Reconstruction by cMEM', 1, size(avg_list, 1)); 
    Results = zeros(size(avg_list, 1) ,size(sHead.Gain, 2), 2, length(sDataIn.Time));


    for iAvg = 1:size(avg_list, 1)
        % Put the data in place
        sDataIn.F   = squeeze(mean(dataMat(avg_list(iAvg,:), : , :), 1));
        sDataIn.Std = []; 
        sDataIn.ChannelFlag = ones(size(sDataIn.F,1), 1);

        % Compute MNE
        sResults = process_nst_cmem('Compute',OPTIONS, ChannelMat, sDataIn);
        % Select HbO and HbR
        sResults = process_nst_cmem('filterResults', sResults, [0, 1, 1, 0]); 
        
        % Store MNE resuls
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
        ResultsMat.DataFile = file_short(AvgFile.FileName);
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
