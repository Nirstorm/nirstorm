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
% Authors: Thomas Vincent, ZhengChen Cai (2015-2016)
%
% TODO: fix baseline as external data input

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
    sProcess.options.TimeSegmentNoise.Type    = 'timewindow';
    sProcess.options.TimeSegmentNoise.Value   = [];
    sProcess.options.TimeSegmentNoise.Class = 'noise_cov';

    sProcess.options.store_sparse_results.Comment = 'Store sparse results';
    sProcess.options.store_sparse_results.Type    = 'checkbox';
    sProcess.options.store_sparse_results.Value   = 0;
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
cortex = in_tess_bst(nirs_head_model.SurfaceFile);
HM.SurfaceFile = nirs_head_model.SurfaceFile;
nb_nodes = size(cortex.Vertices, 1);

%% Load data & baseline
% Load recordings

if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs(1).FileName);
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
    sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
end

% TODO handle raw/data 
nb_samples = length(sDataIn.Time);
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);

if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error(['dMNE source reconstruction works only for dOD data ' ... 
               ' (eg do not use MBLL prior to this process)']);
    return;
end
nb_wavelengths = length(ChannelMat.Nirs.Wavelengths);

measure_tag = 'WL';

% % Default options
% MethodOptions = be_main();
% % Interface to edit options
% MethodOptions = gui_show_dialog('MEM options', @panel_brainentropy, [], [], MethodOptions);
% % Add fields that are not defined by the options of the MEM interface
% if ~isempty(MethodOptions)
%     switch (nirs_head_model.HeadModelType)
%         case {'surface', 'ImageGrid'}
%             MethodOptions.SourceOrient{1} = 'fixed';
%         case 'volume'
%             MethodOptions.SourceOrient{1} = 'free';
%             MethodOptions.flagSourceOrient = [0 0 2 0];
%     end
% end
% % Canceled by user
% if isempty(MethodOptions)
%     return
% end
% % Add options to list
% OPTIONS = process_inverse_2016('Compute');
% OPTIONS.InverseMethod = 'mem';
% OPTIONS = struct_copy_fields(OPTIONS, MethodOptions, 1);
% OPTIONS.DataTypes = {'NIRS'};
% OPTIONS.NoiseCov = [];
OPTIONS.NoiseCov_recompute = sProcess.options.NoiseCov_recompute.Value;
OPTIONS.depth_weigth_MNE = sProcess.options.depth_weightingMNE.Value{1};

%% Run dMNE
pbar = bst_progress('start', 'Reconstruction by wMNE', 'Launching wMNE...');
OPTIONS.Comment = 'wMNE';
OPTIONS.DataFile      = sInputs(1).FileName;
OPTIONS.DataTime      = round(sDataIn.Time,6);
OPTIONS.ResultFile    = [];
OPTIONS.HeadModelFile =  sStudy.HeadModel(sStudy.iHeadModel).FileName;
OPTIONS.FunctionName  = 'wMNE';

if isempty(sProcess.options.TimeSegmentNoise.Value{1})
    OPTIONS.BaselineSegment = [sDataIn.Time(1), sDataIn.Time(end)];
else
    OPTIONS.BaselineSegment  = sProcess.options.TimeSegmentNoise.Value{1};
end   

if isempty(sProcess.options.TimeSegment.Value{1})
    OPTIONS.TimeSegment = [sDataIn.Time(1), sDataIn.Time(end)];
else
    OPTIONS.TimeSegment = sProcess.options.TimeSegment.Value{1};
end    

dOD_sources_wMNE = zeros(nb_nodes, nb_wavelengths, nb_samples);

thresh_dis2cortex = sProcess.options.thresh_dis2cortex.Value{1}.*0.01;
store_sparse_results = sProcess.options.store_sparse_results.Value;

for iwl=1:nb_wavelengths
    swl = [measure_tag num2str(ChannelMat.Nirs.Wavelengths(iwl))];
    selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (sDataIn.ChannelFlag>0)';
    
    OPTIONS.GoodChannel = ones(sum(selected_chans), 1);
    OPTIONS.ChannelFlag   = ones(sum(selected_chans), 1);
    OPTIONS.Channel = ChannelMat.Channel(selected_chans);
    
    gain = get_gains(nirs_head_model.Gain, nirs_head_model.pair_names, iwl, ...
        ChannelMat.Channel, find(selected_chans));
    %% define the reconstruction FOV
    montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel);
    src_locs = montage_info.src_pos;
    det_locs = montage_info.det_pos;
    
    optodes_pos = [src_locs;det_locs];
    % inflate surface 100% to calculate distances to optodes (see BST folder figure_3d.m line 2595)
    iVertices = 1:length(cortex.Vertices);
    % Smoothing factor
    SurfSmoothIterations = ceil(300 * 1 * length(iVertices) / 100000);
    % Calculate smoothed vertices locations
    Vertices_sm = cortex.Vertices;
    Vertices_sm(iVertices,:) = tess_smooth(cortex.Vertices(iVertices,:), 1, SurfSmoothIterations, cortex.VertConn(iVertices,iVertices), 1);
    dis2cortex = pdist2(Vertices_sm,optodes_pos);
    valid_nodes = find(min(dis2cortex,[],2)<thresh_dis2cortex);
 
    
    OPTIONS.cortex_vertices = cortex.Vertices(valid_nodes, :); 
    HM.vertex_connectivity = cortex.VertConn(valid_nodes, valid_nodes);
    
    HM.Gain = gain(:,valid_nodes);
    pbar = bst_progress('text', ['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
             ' nodes that were in VOI and have non-zero sensitivity']);
%     disp(['WL' num2str(iwl) ', kept ' num2str(length(valid_nodes)) ...
%              ' nodes that were in VOI and have non-zero sensitivity']);
    
    OPTIONS.Data = sDataIn.F(selected_chans,:);
    
    
    %% launch dMNE 
    HM.Gain(HM.Gain==0) = min(HM.Gain(HM.Gain>0));
    pbar = bst_progress('text', ['Running wMNE for wavelength #' num2str(iwl) '...']);
    Results = NIRS_wMEM(HM, OPTIONS);
        
    % MNE results
    grid_amp = zeros(nb_nodes, nb_samples);
    sample = be_closest(OPTIONS.TimeSegment([1 end]), OPTIONS.DataTime);
    grid_amp(valid_nodes,sample(1):sample(2)) = Results;
    
    dOD_sources_wMNE(:, iwl, :) = grid_amp;
    
    [sStudy, ResultFile] = add_surf_data(grid_amp, sDataIn.Time, nirs_head_model, ...
                                        ['wMNE sources - ' swl 'nm'], ...
                                        sInputs, sStudy, 'wMNE sources reconstruction', ...
                                        'OD', store_sparse_results);    
    
    OutputFiles{end+1} = ResultFile;
end
pbar = bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');
% Compute dHb
hb_extinctions = nst_get_hb_extinctions(ChannelMat.Nirs.Wavelengths);
% hb_extinctions = [2321.4 1791.7;... % 830HbO 830HbR cm-1.mole-1.L
%     956.89 4932.00];   % 690HbO 690HbR cm-1.mole-1.L
hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L
nirs_hbo_hbr = zeros(nb_nodes, 2, nb_samples);
%nirs_hbo_hbr_mne = zeros(nb_nodes, 2, nb_samples);
for idx=1:length(valid_nodes)
    inode = valid_nodes(idx);
    nirs_hbo_hbr(inode, :, :) = pinv(hb_extinctions) * ...
                                squeeze(dOD_sources_wMNE(inode, :, :));
%     %nirs_hbo_hbr_mne(inode, :, :) = pinv(hb_extinctions) * ...
%                                 squeeze(dOD_sources_MNE(inode, :, :));
end
nirs_hbt = squeeze(sum(nirs_hbo_hbr, 2));

hb_unit_factor = 1e6;
hb_unit = '\mumol.l-1';

hb_types = {'HbO', 'HbR'};
if OPTIONS.depth_weigth_MNE > 0
    function_name='wMNE';
else
    function_name='MNE';
end    
for ihb=1:2
    
    [sStudy, ResultFile] = add_surf_data(squeeze(nirs_hbo_hbr(:,ihb,:)) .* hb_unit_factor,...
                                         sDataIn.Time, nirs_head_model, ...
                                         [function_name ' sources - ' hb_types{ihb}], ...
                                         sInputs, sStudy, [function_name ' sources reconstruction - dHb'], ...
                                         hb_unit, store_sparse_results);    
    OutputFiles{end+1} = ResultFile;
end

[sStudy, ResultFile] = add_surf_data(nirs_hbt .* hb_unit_factor, sDataIn.Time, nirs_head_model, ...
                                     [function_name ' sources - HbT'], ...
                                     sInputs, sStudy, [function_name ' sources reconstruction - dHb'], ...
                                     hb_unit, store_sparse_results);
OutputFiles{end+1} = ResultFile;
pbar = bst_progress('stop', 'Reconstruction by wMNE', 'Finishing...');
% Update Brainstorm database
bst_set('Study', sInputs.iStudy, sStudy);
end

% depth weighted mne

function J = NIRS_wMEM(HM,OPTIONS)
sample_baseline = be_closest(OPTIONS.BaselineSegment([1 end]), OPTIONS.DataTime);
baseline = OPTIONS.Data(:,sample_baseline(1):sample_baseline(2));

%% Normalization by using the mean standard deviation of baseline for each modalities

% Std deviation for every channels on a modality
SD = std(baseline');
% Define the mean standard deviation (MSD) of the present modality
MSD = mean(SD);
%Normalize datas, baseline, gain matrix and EmptyRoom_data with the mean std dev
M = OPTIONS.Data./MSD;
Baseline = baseline./MSD;
G = HM.Gain./MSD;

%% solve wMNE
[n_capt, n_sour] = size(G);
% selection of the data:
sample = be_closest(OPTIONS.TimeSegment([1 end]), OPTIONS.DataTime);
M = M(:,sample(1):sample(2));

if OPTIONS.NoiseCov_recompute 
    Sigma_d    =   inv(diag(diag(real(cov(Baseline')))));
end
param1 = [0.1:0.1:1 1:5:100 100:100:1000]; %the param1 list we tested in wMNE_org

%param1= [ 1/3 ];
% for i=1:size(M,1)
%     param1(i)=1/abs(snr(M(i,:))); 
% end    
pbar = bst_progress('text', 'wMNE, solving MNE by L-curve ... ');

p = OPTIONS.depth_weigth_MNE;
if OPTIONS.NoiseCov_recompute
    Sigma_s = diag(power(diag(G'*Sigma_d*G),p)); 
else
    Sigma_s = diag(power(diag(G'*G),p)); 
end
W = sqrt(Sigma_s);
scale = trace(G*G')./trace(W'*W);       % Scale alpha using trace(G*G')./trace(W'*W)
alpha = param1.*scale;
Fit = [];
Prior = [];
[U,S,V] = svd(G,'econ');
G2 = U*S;
Sigma_s2 = V'*Sigma_s*V;

bst_progress('start', 'wMNE, solving MNE by L-curve ... ' , 'Solving MNE by L-curve ... ', 1, length(param1));
for i = 1:length(param1)
    if OPTIONS.NoiseCov_recompute 
        J = ((G2'*Sigma_d*G2+alpha(i).*Sigma_s2)^-1)*G2'*Sigma_d*M;
    else
        J = ((G2'*G2+alpha(i).*Sigma_s2)^-1)*G2'*M; % Weighted MNE solution
    end
    Fit = [Fit,norm(M-G2*J)];       % Define Fit as a function of alpha
    Prior = [Prior,norm(W*V*J)];          % Define Prior as a function of alpha
    
    bst_progress('inc', 1); 
end
bst_progress('stop');

[~,Index] = min(Fit/max(Fit)+Prior/max(Prior));  % Find the optimal alpha
if OPTIONS.NoiseCov_recompute
    J = ((G'*Sigma_d*G+alpha(Index).*Sigma_s)^-1)*G'*Sigma_d*M;
else
    J = ((G'*G+alpha(Index).*Sigma_s)^-1)*G'*M;    
end

pbar = bst_progress('text', 'wMNE, solving MNE by L-curve ... done');

% figure()
% plot(Prior, Fit,'b.');
% hold on;plot(Prior(Index), Fit(Index),'ro');
% hold off
% xlabel('Norm |WJ|');
% ylabel('Residual |M-GJ|');
% title('L-curve');

end






function gains = get_gains(gain_in, gain_pair_names, iwl, channel_def, selected_chans)

% hash gain pair names
gpair_idx = containers.Map();
for gpic=1:length(gain_pair_names)
    gpair_idx(gain_pair_names{gpic}) = gpic;
end

if ~isempty(strfind(channel_def(1).Name, 'WL'))
    measure_tag = 'WL';
else
    measure_tag = 'Hb';
end

gains = zeros(length(selected_chans), size(gain_in, 3));
for ic=1:length(selected_chans)
    ichan = selected_chans(ic);
    chan_name = channel_def(ichan).Name;
    pair_name = chan_name(1:strfind(chan_name, 'WL')-1);
    if ~isempty(pair_name) && gpair_idx.isKey(pair_name)
        gains(ic, :) = squeeze(gain_in(gpair_idx(pair_name), iwl, :));
    end
end

end

function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
                                              sInputs, sStudy, history_comment, ...
                                              data_unit, store_sparse)
                                          
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
