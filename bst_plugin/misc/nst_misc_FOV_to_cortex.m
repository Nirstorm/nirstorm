function sResults = nst_misc_FOV_to_cortex(sResults, nVertex, valid_vertex, isSaveFactor)
    

    mapping = zeros(nVertex, length(valid_vertex)); 
    for iNode = 1:length(valid_vertex)
        mapping(valid_vertex(iNode), iNode) = 1;
    end
    mapping = sparse(mapping);


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