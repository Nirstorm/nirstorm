function [scalo, WData, OPTIONS, t2] = be_cwt2scalo(WData, OPTIONS, t)

nbF     =   numel(OPTIONS.wavelet.freqs_analyzed);
lim1    =   max(1,be_closest(OPTIONS.wavelet.freqs_analyzed, OPTIONS.ridges.frequency_range(1) )-1);
lim2    =   min(nbF,be_closest(OPTIONS.wavelet.freqs_analyzed, OPTIONS.ridges.frequency_range(2) ) + 1);
            
% Loop on all modalities
scalo   =   {};
Chan    =   numel(WData);
for ii = 1 : Chan
    switch OPTIONS.ridges.method
        case 'inhouse'
            % Compute ridge maps
            [SKLO, www, OPTIONS]    = be_get_scalogram( WData{ii}, OPTIONS);
        case 'LillyOlhede2010'
            [ir,jr]                 =   ridgewalk( full(WData{ii}(lim1:lim2,:))', fliplr(OPTIONS.wavelet.freqs_analyzed(lim1:lim2)) ); 
            SKLO                    =   zeros( lim2-lim1+1,numel(OPTIONS.mandatory.DataTime) );
            idx                     =   sub2ind( size(SKLO), jr(~isnan(ir)), ir(~isnan(ir)) );
            SKLO(idx)               =   1; 
            SKLO([1 end],:)         =   0;
            SCL                     =   zeros( size(WData{ii}) );
            SCL(lim1:lim2,:)        =   SKLO;
    end
    scalo{ii} = sparse(SKLO);   
end

t2 = [];
if ~isempty(t)
    t2 = toc;
    fprintf('\tScalograms computed done in %4.1f seconds\n', t2-t);
end

return