function sResults = nst_misc_FOV_to_cortex(sResults, nVertex, valid_vertex, isSaveFactor)
    
    mapping = sparse(valid_vertex, 1:length(valid_vertex), 1, nVertex, length(valid_vertex));
    
    for iMap = 1:length(sResults)
        if iscell(sResults(iMap).ImageGridAmp)
            sResults(iMap).ImageGridAmp = [ {mapping} sResults(iMap).ImageGridAmp ];
        else
            if isSaveFactor
                sResults(iMap).ImageGridAmp  = {mapping ,  sResults(iMap).ImageGridAmp};
            else
                sResults(iMap).ImageGridAmp  = mapping *  sResults(iMap).ImageGridAmp;
            end
        end
    end

    
end