function sResults_hb = nst_mbll_source(sResults, wavelentghts)
%NST_MBLL_SOURCE Apply the MBLL on the data after source reconstruction 

    assert(length(sResults) == length(wavelentghts), 'Provide the wavelength associated with each map')
    assert(length(sResults) > 1, 'Unable to compute the MNLL with only one wavelentght ')


    
    bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');

    % Compute dHb
    hb_extinctions = nst_get_hb_extinctions(wavelentghts);
    hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L

    if ~iscell(sResults(1).ImageGridAmp)
        dOD_sources =  permute( cat(3, sResults.ImageGridAmp), [1 3 2]);
    else
        isConsistent = 1;
        for iResult = 2:length(sResults)
            isConsistent =  isConsistent && isequal( sResults(1).ImageGridAmp{2}, sResults(iResult).ImageGridAmp{2});
        end

        if isConsistent
            dOD_sources = zeros(size(sResults(1).ImageGridAmp{1}, 1) ,length(sResults),size(sResults(1).ImageGridAmp{1}, 2));
            for iResult = 1:length(sResults)
                dOD_sources(:, iResult, :)  = sResults(iResult).ImageGridAmp{1};
            end
        else
            % If not consistent, go back to full time-course
            dOD_sources = zeros(size(sResults(1).ImageGridAmp{1}, 1) ,length(sResults),size(sResults(1).ImageGridAmp{2}, 2));
            for iResult = 1:length(sResults)
                sResults(iResult).ImageGridAmp = sResults(iResult).ImageGridAmp{1} * sResults(iResult).ImageGridAmp{2};
                dOD_sources(:, iResult, :)  = sResults(iResult).ImageGridAmp;
            end
        end
    end

    Hb_sources = zeros(size(dOD_sources,1), 3, size(dOD_sources,3));
    for inode=1:size(dOD_sources,1)
        Hb_sources(inode, 1:2, :) = pinv(hb_extinctions) * ...
                                    squeeze(dOD_sources(inode, :, :));
    
    end
    Hb_sources(:,3,:) = squeeze(sum(Hb_sources, 2));

    hb_unit_factor = 1e6;
    hb_unit = '\mumol.l-1';
    hb_types = {'HbO', 'HbR','HbT'};

    sResults_hb = repmat(sResults(1), 1, 3);
    for iHb = 1:3

        tmp = strsplit(sResults(end).Comment, '|');

        sResults_hb(iHb).Comment = strjoin( [{strjoin(tmp(1:end-1), '|')} ,   '|' , hb_types{iHb}]);
        sResults_hb(iHb).History = sResults(end).History;
        sResults_hb(iHb).DisplayUnits   = hb_unit;
        sResults_hb(iHb) = bst_history('add', sResults_hb(iHb), 'compute', 'Estimate concentration change');

        if iscell(sResults_hb(iHb).ImageGridAmp )
            sResults_hb(iHb).ImageGridAmp{1} = squeeze(Hb_sources(:,iHb,:)) .* hb_unit_factor;
        else
            sResults_hb(iHb).ImageGridAmp = squeeze(Hb_sources(:,iHb,:)) .* hb_unit_factor;
        end
    end

end

