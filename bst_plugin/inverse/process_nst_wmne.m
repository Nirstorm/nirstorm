function varargout = process_nst_wmne( varargin )

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
% Authors: Edouard Delaire (2022-2025)
% Thomas Vincent, ZhengChen Cai (2015-2016)
%
eval(macro_method);
end

function sProcess = GetDescription() 
    % Description the process
    sProcess.Comment     = '3D reconstruction - wMNE';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', '3D reconstruction'};
    sProcess.Index       = 1501; %0: not shown, >0: defines place in the list of processes
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/NIRSTORM#Inverse_problem_using_MNE';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results'};

    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm',2};
    
    sProcess.options.depth_weightingMNE.Comment = str_pad('Depth weighting factor for MNE',35);
    sProcess.options.depth_weightingMNE.Type    = 'value';
    sProcess.options.depth_weightingMNE.Value   = {0.5, '', 1};
    
    sProcess.options.TimeSegment.Comment = str_pad('Reconstruction Time window:',35);
    sProcess.options.TimeSegment.Type    = 'timewindow';
    sProcess.options.TimeSegment.Value   = [];
    
    
    sProcess.options.NoiseCov_recompute.Comment = 'Compute noise covariance of the baseline MNE';
    sProcess.options.NoiseCov_recompute.Type    = 'checkbox';
    sProcess.options.NoiseCov_recompute.Controller = 'noise_cov';
    sProcess.options.NoiseCov_recompute.Value   = 1;
    

    sProcess.options.TimeSegmentNoise.Comment = str_pad('Baseline Time window:',35);
    sProcess.options.TimeSegmentNoise.Type    = 'baseline';
    sProcess.options.TimeSegmentNoise.Value   = [];
    sProcess.options.TimeSegmentNoise.Class = 'noise_cov';

end

function s = str_pad(s,padsize)
    if (length(s) < padsize)
        s = [repmat('&nbsp;', 1, padsize - length(s)), s];
    end
    s = ['<FONT FACE="monospace">' s '</FONT>'];
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) 

OutputFiles = {};

sStudy = bst_get('Study', sInput.iStudy);

% Load head model
if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider running "NIRS -> Compute head model"');
    return;
end
HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;

% Load channels
ChannelMat = in_bst_channel(sInput.ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error('MNE source reconstruction works only for dOD data (eg do not use MBLL prior to this process)');
    return;
end

%% Load recordings 
sDataIn = in_bst_data(sInput.FileName);

%% Run dMNE
bst_progress('start', 'Reconstruction by MNE', 'Launching MNE...');

OPTIONS  = getOptions(sProcess, HeadModelFile, sInput.FileName);
sResults = Compute(OPTIONS, ChannelMat, sDataIn );

bst_progress('text', 'Saving Results...');
for iMap = 1:length(sResults)
    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' nst_protect_fn_str(sResults(iMap).Comment)]);

    ResultsMat          = sResults(iMap);
    ResultsMat.DataFile = sInput.FileName;
    ResultsMat.Options  = OPTIONS;

    bst_save(ResultFile, ResultsMat, 'v6');
    db_add_data(sInput.iStudy, ResultFile, ResultsMat);

    OutputFiles{end+1} = ResultFile;
end

bst_progress('stop', 'Reconstruction by MNE', 'Finishing...');
end

function OPTIONS = getOptions(sProcess, HeadModelFileName, DataFile)
    sDataIn = in_bst_data(DataFile, {'Time', 'History'});


    OPTIONS.NoiseCov_recompute  = sProcess.options.NoiseCov_recompute.Value;
    OPTIONS.depth_weigth_MNE    = sProcess.options.depth_weightingMNE.Value{1};

    OPTIONS.DataTime      = round(sDataIn.Time,6);
    OPTIONS.History       = sDataIn.History;
    OPTIONS.HeadModelFile = HeadModelFileName;

    if isfield(sProcess.options.TimeSegmentNoise, 'Value') && iscell(sProcess.options.TimeSegmentNoise.Value) && ~isempty(sProcess.options.TimeSegmentNoise.Value) && ~isempty(sProcess.options.TimeSegmentNoise.Value{1})
        OPTIONS.BaselineSegment  = sProcess.options.TimeSegmentNoise.Value{1};
    else    
        OPTIONS.BaselineSegment = [sDataIn.Time(1), sDataIn.Time(end)];
    end   
    
    if isfield(sProcess.options.TimeSegment, 'Value') && iscell(sProcess.options.TimeSegment.Value) && ~isempty(sProcess.options.TimeSegment.Value) && ~isempty(sProcess.options.TimeSegment.Value{1})
        OPTIONS.TimeSegment = sProcess.options.TimeSegment.Value{1};
    else
        OPTIONS.TimeSegment = [sDataIn.Time(1), sDataIn.Time(end)];
    end  

    OPTIONS.thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
end

function sResults = Compute(OPTIONS, ChannelMat, sDataIn )

    nirs_head_model = in_bst_headmodel(OPTIONS.HeadModelFile, 1);
    if ndims(nirs_head_model.Gain) == 3
        nirs_head_model = process_nst_import_head_model('convert_head_model', ChannelMat, nirs_head_model, 0);
    end

    sCortex         = in_tess_bst(nirs_head_model.SurfaceFile);
    nb_nodes        = size(sCortex.Vertices, 1);
    nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);

    %% define the reconstruction FOV
    valid_nodes             = nst_headmodel_get_FOV(ChannelMat, sCortex, OPTIONS.thresh_dis2cortex, sDataIn.ChannelFlag);

    sResults = repmat(db_template('resultsmat'), 1, nb_wavelengths);
    isReconstructed = true(1, nb_wavelengths); 

    for iwl=1:nb_wavelengths
        bst_progress('text', sprintf('Running MNE for wavelength # %d ...', iwl));
        
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
        
        if ~any(selected_chans)
            isReconstructed(iwl) = false;
            continue
        end

        OPTIONS.GoodChannel = ones(sum(selected_chans), 1);
        OPTIONS.ChannelFlag = ones(sum(selected_chans), 1);
        OPTIONS.Channel     = ChannelMat.Channel(selected_chans);
        OPTIONS.Data        = sDataIn.F(selected_chans,:);

        HM.Gain = nirs_head_model.Gain(selected_chans, valid_nodes); 
        HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));

        % MNE results
        [Kernel, tmp] = nst_mne_lcurve(HM, OPTIONS);

        sample      = be_closest(OPTIONS.TimeSegment([1 end]), OPTIONS.DataTime);
        M = zeros(size(OPTIONS.Data));
        M(:,sample(1):sample(2)) = tmp;
        
        sResults(iwl).Comment       = sprintf('MNE sources | %s nm', swl);
        sResults(iwl).ImageGridAmp  = {Kernel, M};
        sResults(iwl).Time          = OPTIONS.DataTime;
        sResults(iwl).Function      = 'MNE';
        sResults(iwl).DisplayUnits  = 'OD';
        sResults(iwl).ChannelFlag   = OPTIONS.ChannelFlag;
        sResults(iwl).HeadModelFile = OPTIONS.HeadModelFile;
        sResults(iwl).HeadModelType = nirs_head_model.HeadModelType;
        sResults(iwl).SurfaceFile   = file_short(nirs_head_model.SurfaceFile);
        sResults(iwl).History       = OPTIONS.History;
        sResults(iwl) = bst_history('add', sResults(iwl), 'compute', 'Compute minimum norm estimate (MNE)');
    end
    
    % Filter reconstructed wavelengh
    sResults = sResults(isReconstructed);

    % Compute MBLL 
    if length(sResults) > 1

        sResults_hb = nst_mbll_source(sResults, ChannelMat.Nirs.Wavelengths(isReconstructed));
        sResults    = [ sResults, sResults_hb];

    end

    isSaveFactor =  1;
    sResults = nst_misc_FOV_to_cortex(sResults, nb_nodes, valid_nodes, isSaveFactor);
end

function [idX] = be_closest(vecGuess, vecRef)
% This function returns the index of the closest value of VECREF to those 
% in VECGUESS

    idX     =   [];
    for ii  =   1 : numel(vecGuess)
        [dum, idX(ii)]  =   min( abs(vecGuess(ii)-vecRef) );     
    end

end