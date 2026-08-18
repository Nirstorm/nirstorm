function sResults_hb = nst_mbll_source(sResults, wavelentghts)
%NST_MBLL_SOURCE Apply the MBLL on the data after source reconstruction 

    assert(length(sResults) == length(wavelentghts), 'Provide the wavelength associated with each map')
    assert(length(sResults) > 1, 'Unable to compute the MNLL with only one wavelentght ')

    bst_progress('text', 'Calculating HbO/HbR/HbT in source space...');

    % Prepare data
    if ~iscell(sResults(1).ImageGridAmp)
        dOD_sources = permute(cat(3, sResults.ImageGridAmp), [3,1,2]);
    else
        isConsistent = 1;
        for iResult = 2:length(sResults)
            isConsistent =  isConsistent && isequal( sResults(1).ImageGridAmp{2}, sResults(iResult).ImageGridAmp{2});
        end

        if isConsistent
            dOD_sources = zeros(length(sResults), size(sResults(1).ImageGridAmp{1}, 1), size(sResults(1).ImageGridAmp{1}, 2));
            for iResult = 1:length(sResults)
                dOD_sources(iResult, :, :)  = sResults(iResult).ImageGridAmp{1};
            end
        % If not consistant, but only 2 maps, we can safely concatenate them.
        elseif length(sResults) == 2 && length(sResults(1).ImageGridAmp) == 2  && length(sResults(2).ImageGridAmp) == 2

            nVertex = size(sResults(1).ImageGridAmp{1}, 1);
            nChannelA = size(sResults(1).ImageGridAmp{1}, 2);
            nChannelB = size(sResults(2).ImageGridAmp{1}, 2);

            sResults(1).ImageGridAmp{1} = [sResults(1).ImageGridAmp{1},  zeros(nVertex, nChannelB)];
            sResults(2).ImageGridAmp{1} = [zeros(nVertex, nChannelA), sResults(2).ImageGridAmp{1}];

            dataA = sResults(1).ImageGridAmp{2};
            dataB = sResults(2).ImageGridAmp{2};
        
            assert( size(dataA, 2) == size(dataB, 2), 'Uncompatible time definition');


            sResults(1).ImageGridAmp{2} = [dataA ;  dataB];
            sResults(2).ImageGridAmp{2} = [dataA ;  dataB];

            dOD_sources = zeros(length(sResults), size(sResults(1).ImageGridAmp{1}, 1), size(sResults(1).ImageGridAmp{1}, 2));
            for iResult = 1:length(sResults)
                dOD_sources(iResult, :, :)  = sResults(iResult).ImageGridAmp{1};
            end

        else  % Otherwise, go back to full time-course

            dOD_sources = zeros(length(sResults), size(sResults(1).ImageGridAmp{1}, 1), size(sResults(1).ImageGridAmp{2}, 2));
            for iResult = 1:length(sResults)
                sResults(iResult).ImageGridAmp = sResults(iResult).ImageGridAmp{1} * sResults(iResult).ImageGridAmp{2};
                dOD_sources(iResult, :, :)  = sResults(iResult).ImageGridAmp;
            end
        end
    end
    
    % Compute MBLLL
    hb_extinctions = nst_get_hb_extinctions(wavelentghts);
    hb_extinctions = hb_extinctions ./10;% mm-1.mole-1.L
    
    Hb_sources = zeros(3, size(dOD_sources,2),  size(dOD_sources,3));
    for inode=1:size(dOD_sources,2)
        Hb_sources(1:2, inode, :) = pinv(hb_extinctions) * ...
                                    squeeze(dOD_sources(1:2, inode, :));
    
    end
    Hb_sources(3, :,:) = squeeze(sum(Hb_sources, 1));

    % Y = pagemtimes(pinv(hb_extinctions), dOD_sources);
    % Hb_sources = [Y; sum(Y, 1)]; 

    % Save output
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
            sResults_hb(iHb).ImageGridAmp{1} = squeeze(Hb_sources(iHb, :, :)) .* hb_unit_factor;
        else
            sResults_hb(iHb).ImageGridAmp = squeeze(Hb_sources(iHb, :, :)) .* hb_unit_factor;
        end
    end

end

