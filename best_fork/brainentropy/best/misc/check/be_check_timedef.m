function [OPTIONS, FLAG] = be_check_timedef(OPTIONS, isRF)

% Get baseline data segment
FLAG    =   0;
OPTIONS.automatic.sampling_rate     =   round( 1 / diff( OPTIONS.mandatory.DataTime([1 2]) ) );

% Fills MSP_data field
if isempty(OPTIONS.optional.TimeSegment)
    OPTIONS.optional.TimeSegment    =   OPTIONS.mandatory.DataTime([1 end]);
end
OPTIONS.optional.TimeSegment        =   be_closest( OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime );
OPTIONS.optional.TimeSegment        =   OPTIONS.mandatory.DataTime(OPTIONS.optional.TimeSegment(1):OPTIONS.optional.TimeSegment(end));

% No baseline
if isempty(OPTIONS.optional.Baseline) && ~any(isRF)
    % Baseline segment definition
    OPTIONS.optional.BaselineTime   =   OPTIONS.mandatory.DataTime;
    OPTIONS.optional.Baseline       =   OPTIONS.mandatory.Data;
    if isempty(OPTIONS.optional.BaselineSegment)
        OPTIONS.optional.BaselineSegment = OPTIONS.mandatory.DataTime(1):1/OPTIONS.automatic.sampling_rate:OPTIONS.mandatory.DataTime(end);
    else
        STb     = be_closest( OPTIONS.optional.BaselineSegment(1), OPTIONS.optional.BaselineTime );
        NDb     = be_closest( OPTIONS.optional.BaselineSegment(end), OPTIONS.optional.BaselineTime );
        OPTIONS.optional.Baseline       =   OPTIONS.optional.Baseline(:, STb:NDb);
        OPTIONS.optional.BaselineTime   =   OPTIONS.optional.BaselineTime(STb:NDb);
    end
    
elseif ~any(isRF)
    % Baseline time definition
    OPTIONS.automatic.BaselineType  =   'independent';
    if isempty(OPTIONS.optional.BaselineTime)
        OPTIONS.optional.BaselineTime       = 1 : size(OPTIONS.optional.Baseline,2);
    end
    
    % Baseline segment definition
    if isempty(OPTIONS.optional.BaselineSegment)
        OPTIONS.optional.BaselineSegment    = OPTIONS.optional.BaselineTime;
        STb     = 1;
        NDb     = size(OPTIONS.optional.Baseline,2);
    else
        STb     = be_closest( OPTIONS.optional.BaselineSegment(1), OPTIONS.optional.BaselineTime );
        NDb     = be_closest( OPTIONS.optional.BaselineSegment(end), OPTIONS.optional.BaselineTime );
    end
    
    % Check
    if STb > NDb || NDb > size(OPTIONS.optional.Baseline,2)
        error('In be_main : bad time definition for the baseline')
    end
    
    % Baseline Channels
    CH = 1 : size(OPTIONS.mandatory.Data,1);
    if isfield(OPTIONS.optional, 'BaselineChannels') && ~isempty(OPTIONS.optional.BaselineChannels)
        CH = [];
        for ii = 1 : size( OPTIONS.mandatory.Data, 1 )
            iD  = strcmp( {OPTIONS.optional.BaselineChannels.Channel.Name}, OPTIONS.optional.Channel(ii).Name );
            if sum(iD); CH(ii) = find(iD); end
        end
    end
    
    if numel(CH) ~= numel(OPTIONS.automatic.GoodChannel)
        error('MEM : Channels of baseline file do not match channels in data file');
    end    
    OPTIONS.optional.Baseline = OPTIONS.optional.Baseline(CH,STb:NDb);
    
end

% Data length check
if strcmp(OPTIONS.mandatory.pipeline, 'cMEM') && strcmp(OPTIONS.clustering.clusters_type, 'static')
    
    % Min data duration : 25samples
    tmSMP   =   be_closest( OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime );
    nSMP    =   diff( tmSMP ) + 1;
    minW    =   OPTIONS.optional.MSP_min_window;
    if nSMP < minW
        % try to expand the window
        neededSMP   =   ceil( (minW - nSMP)/2 );
        
        % available samples : left
        nSMPleft    =   tmSMP(1) - 1;
        
        % available samples : right
        nSMPright   =   numel(OPTIONS.mandatory.DataTime) - tmSMP(end);
        
        if any([nSMPleft nSMPright]<neededSMP)
            FLAG = 1;
            % can't expand the window
            fprintf('\nBEst error:\tdata is too short for stable clustering\n\t\tmust be at least %i samples long\n\n', minW);
        else
            OPTIONS.optional.TimeSegment    =   OPTIONS.mandatory.DataTime(tmSMP(1)-neededSMP : ...
                 tmSMP(end) + neededSMP);
            
            % expand the window
            fprintf('\nBEst warning:\tdata window was too short for stable clustering\n\t\texpanded so it contains %i samples\n', minW);
        end
        
    end
    
end

% Check emptyroom data
if ~isempty( OPTIONS.optional.EmptyRoom_data )
    ERch    =   OPTIONS.optional.EmptyRoom_channels;
    MNch    =   {OPTIONS.optional.Channel.Name};
    if isempty(ERch) && numel(MNch)==size(OPTIONS.optional.EmptyRoom_data,1)
        % same channels as data, assume same order
        OPTIONS.automatic.Emptyroom_data	= OPTIONS.optional.EmptyRoom_data; 
        
    elseif ~isempty(ERch)
        nERD    =   zeros( size(OPTIONS.mandatory.Data,1), size(OPTIONS.optional.EmptyRoom_data,2) );
        nFound  =   [];
        
        % assume only MEG has empty room data
        iMod    =   find(strcmp(OPTIONS.mandatory.ChannelTypes, 'MEG'));
        
        for ii  =   1 : numel(ERch)
            idE     =   find(strcmp(MNch, ERch{ii}));
            if ~isempty(idE)
                nFound  =   [nFound ii]; 
                nERD( idE,: )  =   OPTIONS.optional.EmptyRoom_data(ii,:);
            end
        end
        
        if nFound~=numel(iMod)
            FLAG = 1;
            fprintf('\nBEst error:\tInconsistent channels for emptyroom and data');         
        end
        OPTIONS.automatic.Emptyroom_data    =   nERD;            
               
    else
        FLAG = 1;
        fprintf('\nBEst error:\tNo channels for emptyroom data');
        
    end
end
    
if isfield(OPTIONS.automatic, 'Emptyroom_data')
    OPTIONS.automatic.Emptyroom_time    =   (0:size(OPTIONS.automatic.Emptyroom_data,2)-1)/OPTIONS.automatic.sampling_rate;
end 

return