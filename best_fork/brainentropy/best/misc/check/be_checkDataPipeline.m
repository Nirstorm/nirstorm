function [OPTIONS, isRF]    = be_checkDataPipeline(OPTIONS)

% Check data type
warning('OFF')
isRF = 0;

% Load whole data - not optimized for standalone call
try
    load( be_fullfile(OPTIONS.automatic.iProtocol.STUDIES, OPTIONS.optional.DataFile), 'History', 'F', 'Time');
    OPTIONS.mandatory.DataTime  =   Time;
    OPTIONS.mandatory.Data      =   F;
    isRF    =   strcmp( History(:,2), 'Ridge-based filtering' );

catch
    if ~strcmpi( OPTIONS.optional.DataFile, {'', 'data'} )
        fprintf('In be_main :\tData file cannot be loaded. Check baseline definition\n');
    end
    
end

% Only useful is signal is already ridgefiltered
if ~any(isRF)
    OPTIONS     =   be_PIPELINEoptions( OPTIONS );
    
elseif strcmp(OPTIONS.mandatory.pipeline,'wMEM')
    OPTIONS     =   be_PIPELINEoptions( OPTIONS );
    OPTIONS.mandatory.pipeline    = 'rwMEM';
    
elseif any(isRF)
    OPTIONS.mandatory.pipeline    = 'rMEM';
    OPTIONS.automatic.rMEMfiles   = { strrep(OPTIONS.optional.DataFile, OPTIONS.automatic.iProtocol.STUDIES, '') };
    OPTIONS.mandatory.Data        = [];
    
    % Check if ridges are computed for all modalities 
    iD = find( ~any( ~isnan(F(OPTIONS.automatic.GoodChannel,:)),2 ) );
    if ~isempty(iD)
        NMS = unique( OPTIONS.mandatory.ChannelTypes(iD) );
        NMS = cellfun(@(a) [a ' '], NMS, 'uni', false);
        error(['rMEM : Ridges were not computed for ' NMS{:}]) 
    end
    isRF    =   sum(isRF);
    OPTIONS =   be_PIPELINEoptions( OPTIONS );
    
end

return