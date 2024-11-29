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
    sProcess.OutputTypes = {'data', 'data'};
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

%% Load recordings (TODO handle raw/data) 
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

nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);
measure_tag     = 'WL';
OPTIONS         = getOptions(sProcess,nirs_head_model, sInputs(1).FileName);

%% Run dMNE
bst_progress('start', 'Reconstruction by wMNE', 'Launching wMNE...');
[mapping, dOD_sources_wMNE,Hb_sources] = Compute(OPTIONS,ChannelMat, sDataIn );


bst_progress('text', 'Saving Results...');
if OPTIONS.depth_weigth_MNE > 0
    function_name='wMNE';
else
    function_name='MNE';
end    
for iwl=1:nb_wavelengths
OPTIONS.DataTime  = [];
    swl = [measure_tag num2str(ChannelMat.Nirs.Wavelengths(iwl))];
    [sStudy, ResultFile] = add_surf_data({mapping, squeeze(dOD_sources_wMNE(:,iwl,:))}, sDataIn.Time, nirs_head_model, ...
                                        [function_name ' sources - ' swl 'nm'], ...
                                        sInputs, sStudy, [function_name 'sources reconstruction'], ...
                                        'OD', OPTIONS);    
    
    OutputFiles{end+1} = ResultFile;
end

hb_unit_factor = 1e6;
hb_unit = '\mumol.l-1';
hb_types = {'HbO', 'HbR', 'HbT'};

for ihb=1:3
    
    [sStudy, ResultFile] = add_surf_data({mapping, squeeze(Hb_sources(:,ihb,:)) .* hb_unit_factor},...
                                         sDataIn.Time, nirs_head_model, ...
                                         [function_name ' sources - ' hb_types{ihb}], ...
                                         sInputs, sStudy, [function_name ' sources reconstruction - dHb'], ...
                                         hb_unit, store_sparse_results);    
    OutputFiles{end+1} = ResultFile;
end

bst_progress('stop', 'Reconstruction by wMNE', 'Finishing...');
% Update Brainstorm database
bst_set('Study', sInputs.iStudy, sStudy);
end

function OPTIONS = getOptions(sProcess,HeadModel, DataFile)
    sDataIn = in_bst_data(DataFile);


    OPTIONS.NoiseCov_recompute = sProcess.options.NoiseCov_recompute.Value;
    OPTIONS.depth_weigth_MNE = sProcess.options.depth_weightingMNE.Value{1};

    OPTIONS.Comment = 'wMNE';
    OPTIONS.DataFile      = DataFile;
    OPTIONS.DataTime      = round(sDataIn.Time,6);
    OPTIONS.ResultFile    = [];
    OPTIONS.HeadModelFile =  HeadModel.FileName;
    OPTIONS.FunctionName  = 'wMNE';
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

function [mapping,dOD_sources,Hb_sources] = Compute(OPTIONS,ChannelMat, sDataIn )
    nirs_head_model = in_bst_headmodel(OPTIONS.HeadModelFile);
    cortex = in_tess_bst(nirs_head_model.SurfaceFile);
    
    nb_nodes = size(cortex.Vertices, 1);
    nb_samples = length(sDataIn.Time);
    nb_wavelengths  = length(ChannelMat.Nirs.Wavelengths);

    HM.SurfaceFile = nirs_head_model.SurfaceFile;

    %% define the reconstruction FOV
    thresh_dis2cortex       = OPTIONS.thresh_dis2cortex;
    valid_nodes             = nst_headmodel_get_FOV(ChannelMat, cortex, thresh_dis2cortex, sDataIn.ChannelFlag);

    OPTIONS.MEMpaneloptions.optional.cortex_vertices = cortex.Vertices(valid_nodes, :); 
    HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);

    dOD_sources = zeros(length(valid_nodes), nb_wavelengths, nb_samples);

    for iwl=1:nb_wavelengths
        swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
        selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
        
        OPTIONS.GoodChannel = ones(sum(selected_chans), 1);
        OPTIONS.ChannelFlag = ones(sum(selected_chans), 1);
        OPTIONS.Channel     = ChannelMat.Channel(selected_chans);
        OPTIONS.Data = sDataIn.F(selected_chans,:);

        gain    = nst_headmodel_get_gains(nirs_head_model,iwl, ChannelMat.Channel, find(selected_chans));
        HM.Gain = gain(:,valid_nodes);
        HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));

        bst_progress('text', ['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
                 ' nodes that were in VOI and have non-zero sensitivity']);
    
        %% launch dMNE 
        bst_progress('text', ['Running wMNE for wavelength #' num2str(iwl) '...']);

        % MNE results
        Results = nst_mne_lcurve(HM, OPTIONS);

        grid_amp = zeros(length(valid_nodes), nb_samples);
        sample = be_closest(OPTIONS.TimeSegment([1 end]), OPTIONS.DataTime);
        grid_amp(:,sample(1):sample(2)) = Results;
        
        dOD_sources(:, iwl, :) = grid_amp;
    end

    bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');
    % Compute dHb
    hb_extinctions = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths);
    hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L

    Hb_sources = zeros(length(valid_nodes), 3, nb_samples);
    for inode=1:length(valid_nodes)
        Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                    squeeze(dOD_sources(inode, :, :));
    
    end
    Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));

    mapping = zeros(nb_nodes, length(valid_nodes)); 
    for iNode = 1:length(valid_nodes)
        mapping(valid_nodes(iNode), iNode) = 1;
    end
    mapping = sparse(mapping);
end

function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
                                              sInputs, sStudy, history_comment, ...
                                              data_unit, OPTIONS)
                                          
    if nargin < 8
        data_unit = '';
    end
    
    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                            ['results_NIRS_' protect_fn_str(name)]);
                        
    % ===== CREATE FILE STRUCTURE =====
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    ResultsMat.ImageGridAmp = data;
    ResultsMat.DisplayUnits = data_unit;
    ResultsMat.Time          = time;
    ResultsMat.DataFile      = sInputs.FileName;
    ResultsMat.HeadModelFile = head_model.FileName;
    ResultsMat.HeadModelType = head_model.HeadModelType;
    ResultsMat.Options       = OPTIONS;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = file_short(head_model.SurfaceFile);
    ResultsMat.GridLoc    = [];
    ResultsMat.GridOrient = [];
    ResultsMat.nAvg      = 1;
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


function [idX] = be_closest(vecGuess, vecRef)
% This function returns the index of the closest value of VECREF to those 
% in VECGUESS

idX     =   [];
for ii  =   1 : numel(vecGuess)
    [dum, idX(ii)]  =   min( abs(vecGuess(ii)-vecRef) );     
end

end
