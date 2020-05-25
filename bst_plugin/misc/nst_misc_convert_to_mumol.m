function Y = nst_misc_convert_to_mumol(Y,DisplayUnits)
%NST_MISC_CONVERT_TO_MUMOL Convert Y from DisplayUnits to mumol.l-1
    if strcmp(DisplayUnits, 'mol.l-1')
        Y = Y * 1e6;
    elseif strcmp(DisplayUnits, 'mmol.l-1')
        Y = Y * 1e3;
    elseif strcmp(DisplayUnits, 'mumol.l-1') || strcmp(DisplayUnits, '\mumol.l-1')
        Y = Y * 1;
    else
        if ~isempty(DisplayUnits)
            warning('Cannot interpret data unit: %s.',DisplayUnits);
        else
            warning('Unspecified data unit.');
        end  
    end
end

