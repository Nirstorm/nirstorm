function  prop = nst_get_tissues_optical_properties(tissues,wavelength)

    txt = fileread("tissues_property.json");
    data = jsondecode(txt);
    
    
    % [mua, mus, g, n]
    prop = nan(size(tissues,1), 4 );
    
    for iTissue = 1:size(tissues,1)
    
        tissues_name    = tissues{iTissue,2};
        tissue_idx      = tissues{iTissue,1};
    
        if ~isfield(data,tissues_name)
            error('Unknown tissue type %s', tissues_name);
        end
    
        data_tissue     = data.(tissues_name);
        iWavelength     = find(data_tissue.wavelength == wavelength);
        if isempty(iWavelength)
            error('No optical property for wavelength %d for tissue type %s', wavelength, tissues_name);
        end
    
        prop(1+tissue_idx, 1) = data_tissue.mua(iWavelength);
        prop(1+tissue_idx, 2) = data_tissue.mus(iWavelength);
        prop(1+tissue_idx, 3) = data_tissue.g(iWavelength);
        prop(1+tissue_idx, 4) = data_tissue.n(iWavelength);
    end

end
