function gains = nst_headmodel_get_gains(head_model, unused_iWL, sChannel, selected_channels)

    if ndims(head_model.Gain) == 3
        % Old format. Gain is (nChannel/nWavelength) x nWavelength x nVertex
        all_wavelength      = unique({sChannel.Group});
        selected_wavelength = unique({sChannel(selected_channels).Group});
        
        iWL = find(strcmp(all_wavelength,selected_wavelength));
        assert(unused_iWL == iWL );
        
        gains     = nst_headmodel_get_gains_old(head_model, iWL, sChannel, selected_channels);
    else
        
        gains = head_model.Gain(selected_channels,  :); 
    end

end



function gains = nst_headmodel_get_gains_old(nirs_head_model,iwl, channel_def, selected_chans)

    gain_in         = nirs_head_model.Gain;
    gain_pair_names = nirs_head_model.pair_names;
    
    % hash gain pair names
    gpair_idx = containers.Map();
    for gpic=1:length(gain_pair_names)
        gpair_idx(gain_pair_names{gpic}) = gpic;
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