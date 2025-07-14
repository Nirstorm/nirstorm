function gains = nst_headmodel_get_gains(HeadModel, unused_iWL, sChannel, selected_channels)

    warning('This function is depracted. For more information, consult https://github.com/Nirstorm/nirstorm/pull/266')
    
    % Old head model
    if ndims(HeadModel.Gain) == 3
        HeadModel = process_nst_import_head_model('convert_head_model', sChannel, HeadModel, 0);
    else
        % New Head model. Check that the orientation is applied, or apply
        % it

        if  size(HeadModel.Gain,2) == 3 * length(HeadModel.GridOrient)
            % Apply the fixed orientation to the Gain matrix (normal to the cortex)
            HeadModel.Gain = bst_gain_orient(HeadModel.Gain, HeadModel.GridOrient);
        end
    end
    
    gains = HeadModel.Gain(selected_channels,  :); 

end