function [pair_chan_indexes, pair_sd_ids] = nst_get_pair_indexes_from_names(channel_labels)
% NST_GET_PAIR_INDEXES_FROM_NAMES groups channels by optode pair
%
%   [PAIR_CHAN_INDEXES, PAIR_SD_IDS] = NST_GET_PAIR_INDEXES_FROM_NAMES(CHANNEL_LABELS)
%       CHANNEL_LABELS (cell array of str): 
%           each str is formatted as 'SxDyWLz' or 'SxDyHbt', where:
%               x: source index
%               y: detector index
%               z: wavelength
%               t: Hb type (O, R, T).
%           Consistency is asserted so that:
%               - channels are unique
%               - channel type is homogeneous
%           Warning is issued if:
%               - the number of measures per pair is not homogeneous
%                 eg one pair has only one wavelength
%
%           Examples: S1D2WL685, S01D7WL830, S3D01HbR
%           
%
%        PAIR_CHAN_INDEXES (2D array of int): channel indexes for each pair
%                                             size: (nb_pairs, nb_measures)
%        PAIR_SD_IDS (2D array of int): source and detector ids
%                                       size: (nb_pairs, 2)
%
%   See also NST_UNFORMAT_CHANNELS

[isrcs, idets, measures, channel_type] = nst_unformat_channels(channel_labels);

pair_to_chans = containers.Map();
pair_to_sd = containers.Map();

for ichan=1:length(isrcs)
    src_id = isrcs(ichan);
    det_id = idets(ichan);
    pair_name = sprintf('S%dD%d', src_id, det_id);
    pair_to_sd(pair_name) = [src_id, det_id];
    if pair_to_chans.isKey(pair_name)
        pair_to_chans(pair_name) = [pair_to_chans(pair_name) ichan];
    else
        pair_to_chans(pair_name) = ichan;
    end
end

pair_chan_indexes = cell2mat(pair_to_chans.values');
pair_sd_ids = cell2mat(pair_to_sd.values');
end