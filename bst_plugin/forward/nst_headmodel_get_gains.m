function gains = nst_headmodel_get_gains(nirs_head_model,iwl, channel_def, selected_chans)

gain_in         = nirs_head_model.Gain;
gain_pair_names = nirs_head_model.pair_names;

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