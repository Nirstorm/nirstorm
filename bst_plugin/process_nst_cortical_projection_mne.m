function varargout = process_nst_cortical_projection_mne( varargin )

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
% Authors: Thomas Vincent, Alexis Machado (2018)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Cortical projection - MNE';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'NIRS - wip';
sProcess.Index       = 1206;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'results', 'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 1;


sProcess.options.head_model_fn.Comment = 'Head model file name: ';
sProcess.options.head_model_fn.Type = 'text';
sProcess.options.head_model_fn.Hidden = 1;
sProcess.options.head_model_fn.Value = '';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};
sStudy = bst_get('Study', sInputs.iStudy);

if isempty(sProcess.options.head_model_fn.Value) 
    % Retrieve head model from current study

    if isempty(sStudy.iHeadModel)
        bst_error('No head model found. Consider process "Compute head model from fluence"');
        return;
    end
    head_model_fn = sStudy.HeadModel(sStudy.iHeadModel).FileName;
else
    % Use given head model file
    head_model_fn = sProcess.options.head_model_fn.Value;
end

head_model = in_bst_headmodel(head_model_fn);
ChanneMat = in_bst_channel(sInputs(1).ChannelFile);

if ~strcmp(head_model.HeadModelType, 'surface')
    bst_error('Extraction only works for surface head model');
    return;
end

if ndims(head_model.Gain) ~= 3
    % TODO: better test shape consistency
    bst_error('Bad shape of gain matrix, must be nb_pairs x nb_wavelengths x nb_vertices');
    return;
end

% Load recordings
if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs(1).FileName);
    events = sDataIn.Events;
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
    sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
    events = sDataRaw.F.events;
end

% Separate NIRS channels from others (NIRS_AUX etc.)
[fnirs, fchannel_def, nirs_other, channel_def_other] = ...
    process_nst_mbll('filter_data_by_channel_type', sDataIn.F', ChanneMat, 'NIRS');


% Remove bad channels: they won't enter MBLL computation so no need to keep them
[good_nirs, good_channel_def] = process_nst_mbll('filter_bad_channels', sDataIn.F',...
                                                 ChanneMat, sDataIn.ChannelFlag);
    
[nirs_psig, pair_names, pair_loc, pair_indexes] = process_nst_mbll('group_paired_channels', good_nirs, good_channel_def);
try
    sensitivity_surf = process_nst_import_head_model('get_sensitivity_from_chans', head_model, pair_names);
catch ME  
    if ~isempty(sStudy.iHeadModel)
        bst_warning('Warning', sProcess, sInputs, 'Given head model is not consistent with current montage. Using default one.');
        head_model_fn = sStudy.HeadModel(sStudy.iHeadModel).FileName;
        head_model = in_bst_headmodel(head_model_fn);
        sensitivity_surf = process_nst_import_head_model('get_sensitivity_from_chans', ...
                                                         head_model, pair_names);
    else
        if strcmp(ME.identifier, 'NIRSTORM:HeadmodelMismatch')
            bst_report('Error', sProcess, sInputs, ME.message);
            return;
        else
            rethrow(ME);
        end
    end
end
% nirs_psig (nb_pairs x nb_wavelengths x nb_samples
% nb_nodes = size(sensitivity_surf, 3);
% nb_pairs = size(pair_ichans, 1);

% nb_wavelengths x [HbO, HbR]
ext_coeffs = process_nst_mbll('get_hb_extinctions', ChanneMat.Nirs.Wavelengths);

param.sensors.cov.flag_cov = 1;
param.sensors.cov.window = [sDataIn.Time(1) sDataIn.Time(1)+5];
[pdata_hbo, pdata_hbr] = mfip_inverse_problem(squeeze(nirs_psig(:,1,:)), ...
    squeeze(nirs_psig(:,2,:)), ...
    sDataIn.Time, sensitivity_surf, ...
    ext_coeffs, param);


[sStudy, ResultFile] = add_surf_data(pdata_hbo, sDataIn.Time, ...
    head_model, 'hbo_projected' , ...
    sInputs.iStudy, sStudy,  ...
    'Projected HbO signals');
OutputFiles{end+1} = ResultFile;

[sStudy, ResultFile] = add_surf_data(pdata_hbr, sDataIn.Time, ...
    head_model, 'hbr_projected' , ...
    sInputs.iStudy, sStudy,  ...
    'Projected HbR signals');
OutputFiles{end+1} = ResultFile;

clear pdata_hbo pdata_hbr;

%     OutputFiles{end+1} = ResultFile;
% for iw=1:size(sensitivity_surf, 2)
%     mixing_mat = squeeze(sensitivity_surf(:, iw, :));
%     nz = all(mixing_mat ~= 0, 2);
%     mixing_mat(nz,:) = mixing_mat(nz,:) ./ repmat(sum(mixing_mat(nz,:), 2), 1, nb_pairs);
%     projected_fnirs = mixing_mat * sDataIn.F(pair_ichans(:,iw), :);
%     % TODO: comment name including measure tag
%
%     comment_measure = sprintf('wl%d', ChanneMat.Nirs.Wavelengths(iw));
%
%
%     % label = [sDataIn.Comment '_' comment_measure '_projected'];
%     label = [comment_measure '_projected'];
%     [sStudy, ResultFile] = add_surf_data(projected_fnirs, sDataIn.Time, ...
%         head_model, label , ...
%         sInputs.iStudy, sStudy,  ...
%         'Projected fNIRS signals');
%
%     OutputFiles{end+1} = ResultFile;
%
% end

end

function [pdata_hbo, pdata_hbr] = mfip_inverse_problem(data_w1, data_w2, t, mat_A, ex, param)
% data_wX contains signals (nb_pairs x nb_samples)
% mat_A is lead field matrix (nb_channels x nb_wavelengths x nb_voxels)
% ex contains the extinction coeffs nb_wavelength x [HbO HbR]

% A_wX -> (nb_pairs x nb_voxels)
A_w1 = squeeze(mat_A(:,1,:));
A_w2 = squeeze(mat_A(:,2,:));

nSamples = length(t);

nVox = size(mat_A, 3);

% data_w1=data_w1'; % nSensors by nSamples
% data_w2=data_w2';
data_topo = [data_w1 ; data_w2]; % nb_pairs*2 x nSamples

flag_show_l_curve = 0;
param.inverse.Tikkonov.lambda{1} = 0.5;

g=[A_w1*ex(1,1) A_w1*ex(1,2); A_w2*ex(2,1) A_w2*ex(2,2)];
%g=blkdiag(A_w1,A_w2);
I=speye(size(g,1));
S=svd(g); %[U,S,V] = svd(g);
if flag_show_l_curve & numel(param.inverse.Tikkonov.lambda{1}>1)
    % L curve
    %---------------------
    blocks=1:100:nSamples+1; % +1 because we substract one in the next block of code for the upper limit
    if ~(blocks(end)==nSamples+1)
        blocks(end+1)=nSamples+1;
    end
    
    for iL=1:numel(param.inverse.Tikkonov.lambda{1})
        % M=g'/(g*g'+param.inverse.Tikkonov.lambda{1}(iL)*I);
        M=g'/(g*g'+param.inverse.Tikkonov.lambda{1}(iL)*max(S)*I);
        
        J=zeros(nVox*2,nSamples);
        for iB=1:numel(blocks)-1
            fprintf('block %g\n',iB)
            J(:,blocks(iB):blocks(iB+1)-1)= M*data_topo(:,blocks(iB):blocks(iB+1)-1);
        end
        
        norm_res(iL)=norm(data_topo-g*J,2);
        norm_sol(iL)=norm(J,2);
        
    end
    clear M blocks J; % M cleared??
    
    ismonotonic = ismonotonic(norm_res, 1);
    ismonotonic = ismonotonic(norm_sol, 1);
    
    subplot(1,1,1)
    % plot(log10(norm_res),log10(norm_sol),'o');
    plot(norm_res,norm_sol,'o');
    xlabel('||J||2');
    ylabel('||d-gJ||2');
    axis('tight')
    [k_corner,info] = corner(norm_res,norm_sol);
    hold on
    plot(log10(norm_res(k_corner)),log10(norm_sol(k_corner)),'+');
    param.inverse.Tikkonov.lambda{1}(k_corner)
    
else
    M=g'/(g*g'+param.inverse.Tikkonov.lambda{1}(1)*max(S)*I);
end

%TODO: recommand small amount of vertices
%TODO: use sparse storage?
pdata_hb = M * data_topo;
pdata_hbo = pdata_hb(1:nVox, :);
pdata_hbr = pdata_hb((nVox+1):(2*nVox), :);

clear pdata_hb;
end


function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
    iStudy, sStudy, history_comment)

%% Save a cortical map to brainstorm with given data

ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
    ['results_' protect_fn_str(name)]);

% ===== CREATE FILE STRUCTURE =====
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = name;
ResultsMat.Function      = '';
ResultsMat.ImageGridAmp = data;
ResultsMat.Time          = time;
ResultsMat.DataFile      = [];
if ~isempty(sStudy.iHeadModel)
    ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
end
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
newResult.DataFile      = ''; %sInputs.FileName;
newResult.isLink        = 0;
newResult.HeadModelType = ResultsMat.HeadModelType;
% Add new entry to the database
iResult = length(sStudy.Result) + 1;
sStudy.Result(iResult) = newResult;
% Update Brainstorm database
bst_set('Study', iStudy, sStudy);
end



function [pair_names, pair_loc, pair_ichans, pair_sd_indexes, ...
    src_coords, src_ids, src_ichans, ...
    det_coords, det_ids, det_ichans] = explode_channels(channel_def)
%% Explode channel data according to pairs, sources and detectors
% Args
%    - channel_def: struct
%        Definition of channels as given by brainstorm
%        Used fields: Channel
%
% TOCHECK WARNING: uses containers.Map which is available with matlab > v2008
%
%  Outputs:
%     - pair_names: cell array of str, size: nb_pairs
%         Pair names, format: SXDX
%     - pair_loc: array of double, size: nb_pairs x 3 x 2
%         Pair localization (coordinates of source and detector)
%     - pair_ichans: matrix of double, size: nb_pairs x nb_wavelengths
%         Input channel indexes grouped by pairs
%     - pair_sd_indexes: matrix of double, size: nb_pairs x 2
%         1-based continuours indexes of sources and detectors for each
%         sources.
%     - src_coords:   nb_sources x 3
%         Source coordinates, indexed by 1-based continuous index
%         To access via source ID, as read from pair name:
%             src_coords(src_id2idx(src_ID),:)
%     - src_ids: 1d array of double, size: nb_sources
%         vector of source ids (as used in pair name)
%     - src_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the source belongs (indexed by 1-based
%         continuous index).
%     - det_coords:   nb_detectors x 3
%         Detector coordinates, indexed by 1-based continuous index
%         To access via detector ID, as used in pair name:
%             det_coords(det_id2idx(det_ID),:)
%     - det_ids: 1d array of double, size: max_detector_id (hashing vector)
%         vector of detector ids (as used in pair name)
%     - det_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the detector belongs (indexed by 1-based
%         continuous index).

MT_OD = 1;
MT_HB = 2;

if isfield(channel_def.Nirs, 'Wavelengths')
    nb_measures = length(channel_def.Nirs.Wavelengths);
    measure_type = MT_OD;
else
    nb_measures = length(channel_def.Nirs.Hb);
    measure_type = MT_HB;
end

pair_to_chans = containers.Map();
pair_to_sd = containers.Map();
src_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
src_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ichan=1:length(channel_def.Channel)
    if strcmp(channel_def.Channel(ichan).Type, 'NIRS')
        chan_name = channel_def.Channel(ichan).Name;
        if measure_type == MT_OD
            iwl = strfind(chan_name, 'WL');
            pair_name = chan_name(1:iwl-1);
            wl = str2double(chan_name(iwl+2:end));
            imeasure = channel_def.Nirs.Wavelengths==wl;
        else
            ihb = strfind(chan_name, 'Hb');
            pair_name = chan_name(1:ihb-1);
            imeasure = strcmp(chan_name(ihb:end), channel_def.Nirs.Hb);
        end
        
        if pair_to_chans.isKey(pair_name)
            measures = pair_to_chans(pair_name);
        else
            measures = zeros(1, nb_measures);
        end
        measures(imeasure) = ichan;
        pair_to_chans(pair_name) = measures;
        
        
        [src_id, det_id] = split_pair_name(pair_name);
        pair_to_sd(pair_name) = [src_id, det_id];
        if src_to_chans.isKey(src_id)
            src_to_chans(src_id) = [src_to_chans(src_id) ichan];
        else
            src_to_chans(src_id) = ichan;
            src_coords_map(src_id) = channel_def.Channel(ichan).Loc(:, 1);
        end
        if det_to_chans.isKey(det_id)
            det_to_chans(det_id) = [det_to_chans(det_id) ichan];
        else
            det_to_chans(det_id) = ichan;
            det_coords_map(det_id) = channel_def.Channel(ichan).Loc(:, 2);
        end
        
    end
end

src_coords = cell2mat(src_coords_map.values)';
src_ichans = src_to_chans.values;
src_ids = cell2mat(src_coords_map.keys);

det_coords = cell2mat(det_coords_map.values)';
det_ichans = det_to_chans.values;
det_ids = cell2mat(det_coords_map.keys);

nb_pairs = pair_to_chans.size(1);
pair_names = pair_to_chans.keys;
pair_ichans = zeros(nb_pairs, nb_measures);
pair_loc = zeros(nb_pairs, 3, 2);
pair_sd_indexes = zeros(nb_pairs, 2);
for ipair=1:nb_pairs
    p_indexes = pair_to_chans(pair_names{ipair});
    pair_ichans(ipair, :) = p_indexes;
    pair_loc(ipair, : , :) = channel_def.Channel(pair_ichans(ipair, 1)).Loc;
    sdi = pair_to_sd(pair_names{ipair});
    pair_sd_indexes(ipair, 1) = find(src_ids==sdi(1));
    pair_sd_indexes(ipair, 2) = find(det_ids==sdi(2));
end

end

function [isrc, idet] = split_pair_name(pair_name)
pair_re = 'S([0-9]{1,2})D([0-9]{1,2})';
toks = regexp(pair_name, pair_re , 'tokens');
isrc = str2double(toks{1}{1});
idet = str2double(toks{1}{2});
end


function sfn = protect_fn_str(s)
sfn = strrep(s, ' ', '_');
end
