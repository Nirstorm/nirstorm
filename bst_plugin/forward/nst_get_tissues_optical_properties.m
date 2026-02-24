function  [prop, error_list] = nst_get_tissues_optical_properties(tissues, wavelength, FilesTissuesProperty)
    
    if nargin < 3 || isempty(FilesTissuesProperty) || ~file_exist(FilesTissuesProperty)
        FilesTissuesProperty = fullfile(fileparts(which('nst_get_tissues_optical_properties')),  'tissues_property.json');
    end

    txt  = fileread(FilesTissuesProperty);
    data = jsondecode(txt);
    
    
    % [mua, mus, g, n]
    prop        = nan(size(tissues,1), 4 );
    error_list  = {};

    for iTissue = 1:size(tissues,1)
    
        tissues_name    = tissues{iTissue,2};
        tissue_idx      = tissues{iTissue,1};
    
        if ~isfield(data,tissues_name)
            error_list{end+1} = sprintf('Unknown tissue type %s.', tissues_name);
            continue;
        end
    
        data_tissue     = data.(tissues_name);
        iWavelength     = find(data_tissue.wavelength == wavelength);
        if isempty(iWavelength)
            error_list{end+1} = sprintf('No optical property for wavelength %d for tissue type %s.', wavelength, tissues_name);
            continue;
        end
    
        prop(1+tissue_idx, 1) = data_tissue.mua(iWavelength);
        prop(1+tissue_idx, 2) = data_tissue.mus(iWavelength);
        prop(1+tissue_idx, 3) = data_tissue.g(iWavelength);
        prop(1+tissue_idx, 4) = data_tissue.n(iWavelength);
    end
end
