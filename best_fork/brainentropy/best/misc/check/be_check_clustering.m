function [fault, OPTIONS] = be_check_clustering(OPTIONS, HeadModel)
        
    [nbS, nbT] = size(OPTIONS.optional.clustering.clusters);
    fault = 0;
    
    % format as column vector
    if nbS==1 
        OPTIONS.optional.clustering = OPTIONS.optional.clustering';
        [nbS, nbT] = size(OPTIONS.optional.clustering);
    end
    
    % check nb of sources
    if mod( size(HeadModel.Gain,2), nbS )
        fprintf('MEM error > Faulty manual clustering. Check nb of sources\n\n');
        fault = 1;
    end
    
    % check nb of bins
    DTs = be_closest( OPTIONS.optional.TimeSegment(1), OPTIONS.mandatory.DataTime );
    DTn = be_closest( OPTIONS.optional.TimeSegment(end), OPTIONS.mandatory.DataTime );
    if nbT~=1 && nbT~=DTn-DTs+1
        fprintf('MEM error > Faulty manual clustering. Check nb of time samples\n\n');
        fault = 1;
    end
    
    % check clustering values
    minC = min(OPTIONS.optional.clustering.clusters);
    maxC = max(OPTIONS.optional.clustering.clusters);
    if any(minC<0) 
        fprintf('MEM error > Faulty manual clustering. No negative labels allowed\n\n');
        fault = 1;
    end
    
    for ii = 1 : nbT
        if numel( unique(OPTIONS.optional.clustering.clusters(:,ii)) ) ~= maxC(ii)-minC(ii)+1
            
            % -- Correct discontinuous labeling -- %
            
            parcel          = OPTIONS.optional.clustering.clusters(:,ii);
            [sorted, orig]  = sort(parcel);
            
            % Keep original order
            reso = [1:12294]';
            reso = reso(orig);
            
            % Remove discontinous labels
            left = [0;diff(sorted)];
            left = left - 1;
            left(left<0) = 0;
            left = cumsum(left);
            sorted = sorted - left;
            new(reso) = sorted;
            
            fprintf('MEM warning > Faulty manual clustering CORRECTED. Labeling was discontinuous at sample %i\n\n',ii);
            
        end
    end
    
return