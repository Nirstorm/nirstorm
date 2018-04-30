function [HeadModel, OPTIONS, FLAG] = be_checkio( HeadModel, OPTIONS, verbose )
    
    FLAG = 0;
    
    % Mandatory
    MFs = {'DataTypes', 'ChannelTypes', 'pipeline', 'stand_alone', 'DataTime'};
    if ~any( isfield(OPTIONS.mandatory, MFs) )
        idF = find( ~any( isfield(OPTIONS.mandatory, MFs) ) );
        fprintf('In be_main :\tmandatory field(s) %s undefined.\n', MFs{idF} );
        FLAG = 1;
    end
    if size( OPTIONS.mandatory.ChannelTypes, 1) == 1
        OPTIONS.mandatory.ChannelTypes = OPTIONS.mandatory.ChannelTypes';
    end

    % Check if gain matrix is consistent with data
    if ~isempty(OPTIONS.mandatory.Data)
        % Check data/time match
        if size(OPTIONS.mandatory.DataTime,2) ~= size(OPTIONS.mandatory.Data,2)
            fprintf('In be_main :\tData does not match time definition.\n');
            FLAG = 1;
        end

        % GoodChannel
        if OPTIONS.automatic.stand_alone & ~OPTIONS.automatic.process
            % Check the consistency of # of channels
            GAIN    =   zeros( numel(OPTIONS.mandatory.ChannelTypes), size(HeadModel.Gain(1).matrix,2) );
            nC1     =   numel(OPTIONS.mandatory.ChannelTypes);
            nC2     =   size(OPTIONS.mandatory.Data,1);
            if (nC1~=nC2)
                fprintf('In be_main :\tChannels definition does not match data size.\n');
                FLAG = 1;
            end
            
            % Check the channels of selected modalities
            OPTIONS.automatic.GoodChannel = zeros( nC1,1 );
            for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
                
                % Check if data present for modality ii
                iDc     =   strcmpi(OPTIONS.mandatory.ChannelTypes, OPTIONS.mandatory.DataTypes{ii});
                if ~sum(iDc)
                    fprintf('In be_main :\tNo channels found for data of type %s.\n', OPTIONS.mandatory.DataTypes{ii});
                    FLAG = 1;
                end
                OPTIONS.automatic.GoodChannel = OPTIONS.automatic.GoodChannel + iDc;
                
                % Check if gain present for modality ii
                MODii   =   strcmp(OPTIONS.mandatory.DataTypes{ii}, {HeadModel.Gain.modality});
                if ~any(MODii)
                    fprintf('In be_main :\tNo gain found for data of type %s.\n', OPTIONS.mandatory.DataTypes{ii});
                    FLAG = 1;
                end
                if size(HeadModel.Gain(MODii).matrix,1) ~= sum(iDc)
                    fprintf('In be_main :\tNo match between data and gain channel for modality : %s.\n', OPTIONS.mandatory.DataTypes{ii});
                    FLAG = 1;
                end
                GAIN(iDc,:)   =  HeadModel.Gain(ii).matrix;
                
            end
            OPTIONS.automatic.GoodChannel = find(OPTIONS.automatic.GoodChannel);               
            OPTIONS.ChannelFlag     = ones( nC1,1 );
            HeadModel.Gain          =   GAIN;
            
        end
        
    elseif isempty(OPTIONS.optional.Baseline) && verbose
        fprintf('In be_main :\tNo data nor baseline were found. Aborted.\n');
        FLAG    =   1;
        
    end

    % Check manual clustering
    if ~isempty(OPTIONS.optional.clustering.clusters) 
        [FLAG, OPTIONS] = be_check_clustering(OPTIONS, HeadModel);    
    end
            
return