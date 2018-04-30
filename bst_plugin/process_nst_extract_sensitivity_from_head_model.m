function varargout = process_nst_extract_sensitivity_from_head_model( varargin )

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
% Authors: Thomas Vincent, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract sensitivity surfaces from head model';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1201;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

% Save the new head model
sStudy = bst_get('Study', sInputs.iStudy);

head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
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

[pair_names, pair_loc, pair_ichans, pair_sd_indexes, ...
          src_coords, src_ids, src_ichans, ...
          det_coords, det_ids, det_ichans] = explode_channels(ChanneMat);
nb_sources = size(src_coords, 1);
nb_dets = size(det_coords, 1);
sensitivity_surf = head_model.Gain;
nb_nodes = size(sensitivity_surf, 3);
% Save sensitivities
for iwl=1:size(sensitivity_surf, 2)
    sensitivity_surf_sum = sum(sensitivity_surf(:, iwl, :),  1) ;
    [sStudy, ResultFile] = add_surf_data(repmat(squeeze(sensitivity_surf_sum), [1,2]), [0 1], ...
                                         head_model, ['Summed sensitivities - WL' num2str(iwl)], ...
                                         sInputs.iStudy, sStudy,  ...
                                        'sensitivity imported from MCXlab');
    
    OutputFiles{end+1} = ResultFile;

end

if nb_dets < 100
    time = 1:(nb_sources*100 + nb_dets);
    for iwl=1:size(sensitivity_surf, 2)
        %sens_tmp = zeros(nb_nodes, length(time)) - 1;
        sens_tmp = zeros(nb_nodes, length(time));
        for ipair=1:size(sensitivity_surf, 1)
            [src_id, det_id] = split_pair_name(pair_names{ipair});
            sens_tmp(:, det_id + src_id*100) = squeeze(sensitivity_surf(ipair, iwl, :)); %source id will be minutes, det_it will be seconds
        end
        [sStudy, ResultFile] = add_surf_data(sens_tmp, time, ...
            head_model, ['Sensitivities - WL' num2str(iwl)], ...
            sInputs.iStudy, sStudy, 'Sensitivity import from template'); %TODO better denomitation
        OutputFiles{end+1} = ResultFile;
    end
end

% for iwl=1:size(sensitivity_surf, 2)
%     for ipair=1:size(sensitivity_surf, 1)
%         [sStudy, ResultFile] = add_surf_data(repmat(squeeze(sensitivity_surf(ipair, iwl, :)), [1,2]), [0 1], ...
%                                              newHeadModel, ['Sensitivity - WL' num2str(iwl) ' ' pair_names{ipair}], ...
%                                              iStudy, sStudy, 'pair sensitivity imported from MCXlab');
%         OutputFiles{end+1} = ResultFile;
%     end
% end

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
    ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName; 
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
