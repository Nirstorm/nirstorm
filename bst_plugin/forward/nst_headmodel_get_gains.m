function gains = nst_headmodel_get_gains(head_model, unused_iWL, sChannel, selected_channels)

    warning('This function is depracted. For more information, consult https://github.com/Nirstorm/nirstorm/pull/266')

    if ndims(head_model.Gain) == 3
        head_model = process_nst_import_head_model('convert_head_model', sChannel, head_model, 0);
    end
    
    gains = head_model.Gain(selected_channels,  :); 
end