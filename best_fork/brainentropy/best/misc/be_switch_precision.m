function [OPTIONS] = be_switch_precision( OPTIONS, func )

myhndl = eval(['@' func]);

nMod = numel(OPTIONS.automatic.Modality);

for ii =  1: nMod
    OPTIONS.automatic.Modality(ii).data         =   myhndl( OPTIONS.automatic.Modality(ii).data );
    OPTIONS.automatic.Modality(ii).baseline     =   myhndl( OPTIONS.automatic.Modality(ii).baseline );
    OPTIONS.automatic.Modality(ii).covariance   =   myhndl( OPTIONS.automatic.Modality(ii).covariance );
    
    if isfield(OPTIONS.automatic.Modality(ii), 'idata')
        OPTIONS.automatic.Modality(ii).idata    =   myhndl( OPTIONS.automatic.Modality(ii).idata );
    end
    
    if ~isempty( fieldnames(OPTIONS.automatic.Modality(ii).mspDATA) ) 
        OPTIONS.automatic.Modality(ii).mspDATA.F=   myhndl( OPTIONS.automatic.Modality(ii).mspDATA.F );
    end
    
    if ~isempty( fieldnames(OPTIONS.automatic.Modality(ii).gain_struct) )
        OPTIONS.automatic.Modality(ii).gain_struct.Gn       =   myhndl( OPTIONS.automatic.Modality(ii).gain_struct.Gn );
        OPTIONS.automatic.Modality(ii).gain_struct.U        =   myhndl( OPTIONS.automatic.Modality(ii).gain_struct.U );
        OPTIONS.automatic.Modality(ii).gain_struct.lambda   =   myhndl( OPTIONS.automatic.Modality(ii).gain_struct.lambda );
    end
    
end

return