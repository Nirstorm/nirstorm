function montage_info = nst_montage_info_from_bst_channels(channel_def, channel_flags)
%% Explode channel data according to pairs, sources and detectors
% Args
%    - channel_def: struct array
%        Definition of channels as given by brainstorm (see db_template('channeldesc'))
%   [- channel_flags]: 1D array of int
%        Allow channel filtering (only keep channel_flags==1)
% TOCHECK WARNING: uses containers.Map which is available only from matlab > v2008
%
%  Outputs: struct with following fields: 
%     - pair_names: cell array of str, size: nb_pairs
%         Pair names, format: SXDX
%     - pair_loc: array of double, size: nb_pairs x 3 x 2
%         Pair localization (coordinates of source and detector)
%     - pair_ichans: matrix of double, size: nb_pairs x nb_wavelengths
%         Input channel indexes grouped by pairs
%     - pair_sd_indexes: matrix of double, size: nb_pairs x 2
%         1-based continuous indexes of sources and detectors for each
%         pair.
%     - src_pos:   nb_sources x 3
%         Source coordinates, indexed by 1-based continuous index
%         To access via source ID use field src_ids as:
%             src_pos(src_ids == src_ID,:)
%     - src_ids: 1d array of double, size: nb_sources
%         vector of source ids (as used in pair name)
%     - src_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the source belongs (indexed by 1-based
%         continuous index).
%     - det_pos:   nb_detectors x 3
%         Detector coordinates, indexed by 1-based continuous index
%         To access via detector ID, use field det_ids as:
%             det_pos(det_ids == det_ID,:)
%     - det_ids: 1d array of double, size: max_detector_id (hashing vector)
%         vector of detector ids (as used in pair name)
%     - det_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the detector belongs (indexed by 1-based
%         continuous index).
%
% TODO: test channel flags

if nargin < 2
    channel_flags = ones(1, length(channel_def));
end

[isrcs, idets, chan_measures, measure_type] = nst_unformat_channels({channel_def.Name});
measure_types = nst_measure_types();
all_measures = unique(chan_measures(~isnan(chan_measures)));
nb_measures = length(all_measures);

%TODO: reuse already parsed isrcs / idets / measures

pair_to_chans = containers.Map();
pair_to_sd = containers.Map();
src_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
src_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ichan=1:length(channel_def)
    if strcmp(channel_def(ichan).Type, 'NIRS') && channel_flags(ichan)==1
        chan_name = channel_def(ichan).Name;
        switch measure_type
            case measure_types.WAVELENGTH
                iwl = strfind(chan_name, 'WL');
                pair_name = chan_name(1:iwl-1);
                wl = str2double(chan_name(iwl+2:end));
                imeasure = all_measures==wl;
            case measure_types.HB
                ihb = strfind(chan_name, 'Hb');
                pair_name = chan_name(1:ihb-1);
                imeasure = strcmp(chan_name(ihb:end), all_measures);
            otherwise
                error('Bad measure type: %d', measure_type);
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
            src_coords_map(src_id) = channel_def(ichan).Loc(:, 1);
        end
        if det_to_chans.isKey(det_id)
            det_to_chans(det_id) = [det_to_chans(det_id) ichan];
        else
            det_to_chans(det_id) = ichan;
            det_coords_map(det_id) = channel_def(ichan).Loc(:, 2);
        end
    
    end
end

montage_info.src_pos = cell2mat(src_coords_map.values)';
montage_info.src_ichans = src_to_chans.values;
montage_info.src_ids = cell2mat(src_coords_map.keys);

montage_info.det_pos = cell2mat(det_coords_map.values)';
montage_info.det_ichans = det_to_chans.values;
montage_info.det_ids = cell2mat(det_coords_map.keys);

nb_pairs = pair_to_chans.size(1);
pair_names = pair_to_chans.keys;
pair_ichans = zeros(nb_pairs, nb_measures);
pair_loc = zeros(nb_pairs, 3, 2);
pair_sd_indexes = zeros(nb_pairs, 2);
for ipair=1:nb_pairs
    p_indexes = pair_to_chans(pair_names{ipair});
    pair_ichans(ipair, :) = p_indexes;
    pair_loc(ipair, : , :) = channel_def(pair_ichans(ipair, 1)).Loc;
    sdi = pair_to_sd(pair_names{ipair});
    pair_sd_indexes(ipair, 1) = find(montage_info.src_ids==sdi(1));
    pair_sd_indexes(ipair, 2) = find(montage_info.det_ids==sdi(2));
end 

montage_info.pair_sd_indexes = pair_sd_indexes;
montage_info.pair_loc = pair_loc;
montage_info.pair_ichans = pair_ichans;
montage_info.pair_names = pair_names;
end

function [isrc, idet] = split_pair_name(pair_name)
formats = nst_get_formats();
toks = regexp(pair_name, formats.pair_re , 'names');
isrc = str2double(toks.src_id);
idet = str2double(toks.det_id);
end
