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
% Authors: Edouard Delaire (2022-2024)
% Thomas Vincent, ZhengChen Cai (2015-2016)
%
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Source reconstruction - wMNE';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1501; %0: not shown, >0: defines place in the list of processes
    sProcess.isSeparator = 0;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results', 'results'};

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
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

sStudy = bst_get('Study', sInputs.iStudy);
%% Load head model
if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider running "NIRS -> Compute head model"');
    return;
end

nirs_head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
nirs_head_model.FileName = sStudy.HeadModel(sStudy.iHeadModel).FileName;

%% Load recordings 

if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs(1).FileName);
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
end

ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error(['dMNE source reconstruction works only for dOD data ' ... 
               ' (eg do not use MBLL prior to this process)']);
    return;
end



%% Run dMNE
bst_progress('start', 'Reconstruction by MNE', 'Launching MNE...');

OPTIONS  = getOptions(sProcess,nirs_head_model, sInputs(1).FileName);
sResults = Compute(OPTIONS,ChannelMat, sDataIn );


bst_progress('text', 'Saving Results...');

for iMap = 1:length(sResults)

    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName),  ['results_NIRS_' protect_fn_str(sResults(iMap).Comment)]);

    ResultsMat = sResults(iMap);
    ResultsMat.DataFile   = sInputs.FileName;
    ResultsMat.HeadModelFile = OPTIONS.HeadModelFile;
    ResultsMat.Options       = OPTIONS;
    ResultsMat.SurfaceFile   = file_short(nirs_head_model.SurfaceFile);

    bst_save(ResultFile, ResultsMat, 'v6');
    db_add_data( sInputs.iStudy, ResultFile, ResultsMat);

    OutputFiles{end+1} = ResultFile;

end

bst_progress('stop', 'Reconstruction by MNE', 'Finishing...');

end

function OPTIONS = getOptions(sProcess,HeadModel, DataFile)
    sDataIn = in_bst_data(DataFile);


    OPTIONS.NoiseCov_recompute = sProcess.options.NoiseCov_recompute.Value;
    OPTIONS.depth_weigth_MNE = sProcess.options.depth_weightingMNE.Value{1};

    OPTIONS.Comment       = 'MNE';
    OPTIONS.DataFile      = DataFile;
    OPTIONS.DataTime      = round(sDataIn.Time,6);
    OPTIONS.ResultFile    = [];
    OPTIONS.HeadModelFile =  file_short(HeadModel.FileName);
    OPTIONS.FunctionName  = 'MNE';

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

function sResults = Compute(OPTIONS,ChannelMat, sDataIn )

    nirs_head_model = in_bst_headmodel(OPTIONS.HeadModelFile);
    sCortex         = in_tess_bst(nirs_head_model.SurfaceFile);
    
    nb_nodes = size(sCortex.Vertices, 1);
    nb_samples = length(sDataIn.Time);
    nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);


    %% define the reconstruction FOV
    valid_nodes             = nst_headmodel_get_FOV(ChannelMat, sCortex, OPTIONS.thresh_dis2cortex, sDataIn.ChannelFlag);
    HM.SurfaceFile = nirs_head_model.SurfaceFile;
    HM.vertex_connectivity  = sCortex.VertConn(valid_nodes, valid_nodes);
    OPTIONS.MEMpaneloptions.optional.cortex_vertices = sCortex.Vertices(valid_nodes, :); 



    sResults = repmat(db_template('resultsmat'), 1, nb_wavelengths);
    isReconstructed = true(1, nb_wavelengths); 

    for iwl=1:nb_wavelengths
        
        bst_progress('text', sprintf('Running wMNE for wavelength # %d ...', iwl));
        
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

        gain    = nst_headmodel_get_gains(nirs_head_model, iwl, ChannelMat.Channel, find(selected_chans));
        HM.Gain = gain(:,valid_nodes);
        HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));

        % MNE results
        Results = nst_mne_lcurve(HM, OPTIONS);

        grid_amp    = zeros(length(valid_nodes), nb_samples);
        sample      = be_closest(OPTIONS.TimeSegment([1 end]), OPTIONS.DataTime);
        grid_amp(:,sample(1):sample(2)) = Results;
        


        sResults(iwl).Comment       = sprintf('MNE sources | %s nm', swl);
        sResults(iwl).ChannelFlag   = OPTIONS.ChannelFlag;
        sResults(iwl).ImageGridAmp  = grid_amp;
        sResults(iwl).Time          = OPTIONS.DataTime;
        sResults(iwl).DisplayUnits  = 'OD';
        sResults(iwl).HeadModelType = file_short(nirs_head_model.HeadModelType);
        sResults(iwl).Function      = 'MNE';

    end
    % Filter reconstructed wavelengh
    sResults = sResults(isReconstructed);

    % Compute MBLL 
    if length(sResults) > 1

        sResults_hb = nst_mbll_source(sResults, ChannelMat.Nirs.Wavelengths(isReconstructed));
        sResults    = [ sResults, sResults_hb];

    end

    mapping = zeros(nb_nodes, length(valid_nodes)); 
    for iNode = 1:length(valid_nodes)
        mapping(valid_nodes(iNode), iNode) = 1;
    end

    mapping = sparse(mapping);
    for iMap = 1:length(sResults)
        sResults(iMap).ImageGridAmp  = {mapping ,  sResults(iMap).ImageGridAmp};

    end
    
end

function [idX] = be_closest(vecGuess, vecRef)
% This function returns the index of the closest value of VECREF to those 
% in VECGUESS

idX     =   [];
for ii  =   1 : numel(vecGuess)
    [dum, idX(ii)]  =   min( abs(vecGuess(ii)-vecRef) );     
end

end

function sfn = protect_fn_str(s)
    sfn = strrep(s, ' | ', '--');
    sfn = strrep(s, ' : ', '--');
    sfn = strrep(s, ' :', '--');
    sfn = strrep(s, ' ', '_');
end