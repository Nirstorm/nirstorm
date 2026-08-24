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
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    
    sProcess.options.label0.Comment = 'Bootstraps options:';
    sProcess.options.label0.Type = 'label';

    
    sProcess.options.SNR.Comment = str_pad('SNR range', 35);
    sProcess.options.SNR.Type    = 'baseline';
    sProcess.options.SNR.Value   =  [];
    sProcess.options.SNR.Group   = 'MNE';


    sProcess.options.replacement.Comment = 'Do bootstrap with replacement ?';
    sProcess.options.replacement.Type    = 'checkbox';
    sProcess.options.replacement.Value   = 1;

    sProcess.options.combination.Comment = 'Combination value #:';
    sProcess.options.combination.Type    = 'value';
    sProcess.options.combination.Value   = {16, '', 0};
    
    sProcess.options.data_save.Comment = 'Data save comment';
    sProcess.options.data_save.Type    = 'text';
    sProcess.options.data_save.Value = '';

    sProcess.options.label1.Comment = 'MNE options:';
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
% ===== GET OPTIONS =====
% Get options values
OutputFiles = {};

sStudy = bst_get('Study', sInputs.iStudy);
HeadModelFileName = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
HeadModelFileName.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);


for iFile = 1:(length(sInputs))
    DataMat_dOD{iFile} = in_bst_data(sInputs(iFile).FileName, 'F', 'ChannelFlag', 'History', 'Time');
end

ChannelFlag = DataMat_dOD{1,1}.ChannelFlag & strcmp({ChannelMat.Channel.Type},'NIRS')'; 

%baselinedur = [sProcess.options.TimeSegmentNoise.Value{1}, sProcess.options.TimeSegmentNoise.Value{2}];
%SNRdur = [sProcess.options.SNR.Value{1}, sProcess.options.SNR.Value{2}];
Time = round(DataMat_dOD{1,1}.Time,1);

ibaseline = panel_time('GetTimeIndices', DataMat_dOD{1,1}.Time, baselinedur);
iSNR = panel_time('GetTimeIndices', DataMat_dOD{1,1}.Time, SNRdur);

if length(sInputs) <= 20  &&  ~sProcess.options.replacement.Value
    avg_list = nchoosek([1:(length(sInputs))],sProcess.options.combination.Value{1});%Cnk
    n_perm = size(avg_list,1);
else
    n_perm = 5000;
    avg_list = zeros(n_perm,sProcess.options.combination.Value{1});
end

bst_progress('start', 'Bootstraping', 'Computing all combinations', 0, length(avg_list)); 

for iAvg = 1:n_perm
    if length(sInputs) > 20 || sProcess.options.replacement.Value
        avg_list(iAvg,:) =  randsample(length(sInputs),sProcess.options.combination.Value{1}, sProcess.options.replacement.Value ); 
        % avg_list(iavg,:) = randperm(length(sInputs),sProcess.options.combination.Value{1}); 
    end
    trial = cellfun(@(c) c.F,DataMat_dOD(avg_list(iAvg,:)),'UniformOutput',false);

    trials = mean(cat(3,trial{:}),3);
    trials = trials(ChannelFlag==1,:);

    ave_trials_690 = mean(trials(1:2:end,:),1);
    ave_trials_830 = mean(trials(2:2:end,:),1);
    trials_690 = trials(1:2:end,:);
    trials_830 = trials(2:2:end,:);

    trials_std_690 = std(trials_690(:,ibaseline),[],2);
    trials_std_830 = std(trials_830(:,ibaseline),[],2);
    
    SNR_690(iAvg) = max(max(abs(trials_690(:, iSNR)))./mean(trials_std_690));
    SNR_830(iAvg) = max(max(abs(trials_830(:,iSNR)))./mean(trials_std_830));
    avg_SNR_690(iAvg) = abs(min(ave_trials_690(iSNR)))./std(ave_trials_690(iSNR));
    avg_SNR_830(iAvg) = abs(max(ave_trials_830(iSNR)))./std(ave_trials_830(iSNR));
    
    bst_progress('inc', 1);

end

%% generate the avg list of each wavelength around the median amp (e.g.101 resamples)

nAverage = min(50,floor(length(avg_list)/2-1));


% SNR_all = (SNR_690+SNR_830)./2;
% SNR_all = min(SNR_690,SNR_830);
quantiles690 = quantile(SNR_690',[0.5 0.6 0.7 0.8 0.9]);
quantiles830 = quantile(SNR_830',[0.5 0.6 0.7 0.8 0.9]);

%this is targeting median SNR for both wavelength
SNR_all  = sqrt((SNR_690-  quantiles690(1)).^2 + (SNR_830 - quantiles830(1)).^2); 

[SNR_sorted,SNR_all_list] = sort(SNR_all);

%target_SNR = median(SNR_sorted);
target_SNR = min(SNR_sorted);

%target_SNR = max( median(SNR_690),median(SNR_830));

[~,idx_list] =  min( abs(SNR_sorted - target_SNR) );

%avg_all_list = SNR_all_list(idx_list-nAverage:min(length(avg_list), min(length(SNR_all), idx_list+nAverage)));
avg_all_list = SNR_all_list(idx_list:min(length(avg_list), min(length(SNR_all), idx_list+2*nAverage)));


SNR_690_s   = SNR_690(avg_all_list);
SNR_830_s   = SNR_830(avg_all_list);
SNR_s       = (SNR_690_s+SNR_830_s)/2;

figures_process_MNE(avg_list,avg_all_list,sProcess,sInputs,SNR_690,SNR_690_s,SNR_830,SNR_830_s,SNR_all,SNR_s,quantiles830,quantiles690);

% Genrate median average
sFiles = bst_process('CallProcess', 'process_average', sInputs(avg_list(SNR_all_list(idx_list),:)), [], ...
    'avgtype',       5, ...  % By trial group (folder average)
    'avg_func',      7, ...  % Arithmetic average + Standard error
    'weighted',      0, ...
    'keepevents',    1);

% Process: Add tag:  
AvgFile = bst_process('CallProcess', 'process_add_tag', sFiles, [], ...
    'tag',           sprintf( '| median SNR = %.2f',median(SNR_690_s)) , ...
    'output',        1);  % Add to file name


sStudy = bst_get('Study', sInputs.iStudy);
if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider running "NIRS -> Compute head model"');
    return;
end
HeadModelFileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;
sHead = in_bst_headmodel(HeadModelFileName, 1);

sDataIn = in_bst_data(AvgFile.FileName);
ChannelMat = in_bst_channel(AvgFile.ChannelFile);
OPTIONS = process_nst_wmne('getOptions', sProcess, HeadModelFileName, AvgFile.FileName);

bst_progress('start', 'Reconstruction by cMEM', 'Launching cMEM...');


Results = zeros(length(avg_all_list) ,size(sHead.Gain, 2), 2, length(sDataIn.Time));

bst_progress('start', 'Bootstraping', 'Computing all combinations', 0, length(avg_all_list)); 
for iAvg = 1:max(length(avg_all_list), 10)
    
    trial = cellfun(@(c) c.F, DataMat_dOD(avg_list(avg_all_list(iAvg),:)),'UniformOutput',false);

    trials = mean(cat(3, trial{:}),3);
    trials = trials(ChannelFlag==1,:);

    sDataIn.F = trials;

    [sResults] = process_nst_wmne('Compute', OPTIONS, ChannelMat, sDataIn );
    sResults = process_nst_wmne('filterResults', sResults, [0, 1, 1, 0]);

    Results(iAvg,:,1,:) = bst_multiply_cellmat(sResults(1).ImageGridAmp);
    Results(iAvg,:,2,:) = bst_multiply_cellmat(sResults(2).ImageGridAmp);
    
    bst_progress('inc', 1);
end 

ResultsAvg = squeeze(mean(Results));
ResultsSD  = squeeze(std(Results));

hb_unit_factor = 1e6;
hb_unit = '\mumol.l-1';
hb_types = {'HbO', 'HbR'};
for ihb = 1:size(ResultsAvg, 2)
    [sStudy, ResultFile] = add_surf_data(squeeze(ResultsAvg(:,ihb,:)) .* hb_unit_factor,...
                                         squeeze(ResultsSD(:,ihb,:)) .* hb_unit_factor,...
                                         sDataIn.Time, HeadModelFileName, ...
                                         sprintf('MNE sources (bootstrap) | %s', hb_types{ihb}), ...
                                         AvgFile, sStudy, '', ...
                                         length(avg_all_list),...
                                         hb_unit, 0);    
    OutputFiles{end+1} = ResultFile;
end

end

function [sStudy, ResultFile] = add_surf_data(data,Std, time, head_model, name, ...
                                              sInputs, sStudy, history_comment, ...
                                              nAvg, data_unit, store_sparse)
                                          
    if nargin < 8
        data_unit = '';
    end
    
    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                            ['results_NIRS_' protect_fn_str(name)]);
                        
    % ===== CREATE FILE STRUCTURE =====
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    if store_sparse
        ResultsMat.ImageGridAmp = sparse(data); %TODO TOCHECK with FT: sparse data seem not well handled. Eg while viewing (could not reproduce)
    else
        ResultsMat.ImageGridAmp = data;
    end
    ResultsMat.Std = Std;
    ResultsMat.DisplayUnits = data_unit;
    ResultsMat.Time          = time;
    ResultsMat.DataFile      = sInputs.FileName;
    ResultsMat.HeadModelFile = head_model.FileName;
    ResultsMat.HeadModelType = head_model.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = file_short(head_model.SurfaceFile);
    ResultsMat.GridLoc    = [];
    ResultsMat.GridOrient = [];
    ResultsMat.nAvg      = nAvg;
    % History
    ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
    % Save new file structure
    bst_save(ResultFile, ResultsMat, 'v6');
    % ===== REGISTER NEW FILE =====
    % Create new results structure
    newResult = db_template('results');
    newResult.Comment       = name;
    newResult.FileName      = file_short(ResultFile);
    newResult.DataFile      = sInputs.FileName;
    newResult.isLink        = 0;
    newResult.HeadModelType = ResultsMat.HeadModelType;
    % Add new entry to the database
    iResult = length(sStudy.Result) + 1;
    sStudy.Result(iResult) = newResult;
    % Update Brainstorm database
    bst_set('Study', sInputs.iStudy, sStudy);
end



function sfn = protect_fn_str(s)
sfn = strrep(s, ' ', '_');
end

function figures_process_MNE(avg_list,avg_all_list,sProcess,sInputs,SNR_690,SNR_690_s,SNR_830,SNR_830_s,SNR_all,SNR_s,quantiles830,quantiles690)
% TODO: replace hist with histogram
figure
subplot(2,2,1);
histogram(avg_list(avg_all_list,:),'BinMethod','integers')
line(xlim(gca), length(avg_all_list)*(sProcess.options.combination.Value{1}/length(sInputs))*[1,1],'Color','red','LineStyle','--')
xlabel('Trials');
title('Occurance of each trial in the final result')

subplot(2,2,2);
hist(SNR_690,size(SNR_690,2)./10);
line(median(SNR_690_s)*[1,1],ylim(gca),'Color','red','LineStyle','--')
line(min(SNR_690_s)*[1,1],ylim(gca),'Color','blue','LineStyle','--')
line(max(SNR_690_s)*[1,1],ylim(gca),'Color','blue','LineStyle','--')
xlim([0 43])
title('SNR for \lambda = 690');

subplot(2,2,3);
hist(SNR_830,size(SNR_830,2)./10);
line(median(SNR_830_s)*[1,1],ylim(gca),'Color','red','LineStyle','--')
line(min(SNR_830_s)*[1,1],ylim(gca),'Color','blue','LineStyle','--')
line(max(SNR_830_s)*[1,1],ylim(gca),'Color','blue','LineStyle','--')
xlim([0 40])
title('SNR for \lambda = 830');

subplot(2,2,4);
hist((SNR_690+SNR_830)./2,size(SNR_all,2)./10);
line(median(SNR_s)*[1,1],ylim(gca),'Color','red','LineStyle','--')
line(min(SNR_s)*[1,1],ylim(gca),'Color','blue','LineStyle','--')
line(max(SNR_s)*[1,1],ylim(gca),'Color','blue','LineStyle','--')
xlim([0 40])
title('Average SNR for \lambda = {690, 830}');

figure;
scatterhist(SNR_690,SNR_830,'Location','NorthEast', 'Direction','out'); hold on;
scatter(SNR_690_s, SNR_830_s,'filled','r')
yline(quantiles830 ); xline(quantiles690)
sgtitle('SNR for \lambda = {690, 830}')
end