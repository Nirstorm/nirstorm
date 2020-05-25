function model = nst_glm_add_regressors(model,Regressor_type,varargin)
% Call exemple 
% model = add_regressors(Model, "event", Events, basis_choice, hrf_duration)
% model = add_regressors(Model, "external_input", external_input)
% model = add_regressors(Model, "channel", channel_names, [convolve_with_hrf=0, basis_choice, hrf_duration]) Note: if convolve_with_hrf=1, then convolute the regressor with HRF like for event based regressor. Typically if the channel is AUX that contain the paradigm of the experiment
% model = add_regressors(Model, "constant") Add a regressor that just contains 1.
% model = add_regressors(Model, "DCT", frequences, names)


    switch Regressor_type
        case'event'
            
            events=varargin{1};
            if nargin < 4 
                hrf_types   = process_nst_glm_fit('get_hrf_types');
                basis_choice= hrf_types.CANONICAL;
            else
                basis_choice=varargin{2};
            end
            
            if nargin < 5 
                hrf_duration= 25;
            else
                hrf_duration=varargin{3};
            end
            model=nst_glm_add_event_regressors(model,events,basis_choice,hrf_duration);
        case 'external_input'
            sInput_ext=varargin{1};
            if nargin < 4
                hb_types = {'HbO', 'HbR', 'HbT'};
            else
                hb_types=varargin{2};
            end    
            
            model=nst_glm_add_ext_input_regressors(model,sInput_ext,hb_types); 
        case 'channel'
            
            
        case 'constant'
            model=nst_glm_add_constant_regressors(model);
        
        case 'DCT'
            frequences=varargin{1};
            names=varargin{2};
            model=nst_glm_add_DCT_regressors(model,frequences,names);
    end    


end

function model=nst_glm_add_event_regressors(model,events,hrf_type, hrf_duration)
    % Todo, maybe make sure that event are at the egening of the design
    % matrix. Can be forced by changing model.X= [model.X X_event] to
    % model.X= [X_event model.X]; 
    
    % reconstruct time course
    time=model.time;
    dt = time(2)-time(1);
    hrf_time = 0:dt:hrf_duration; 
    
    if hrf_time(end) ~= hrf_duration
        warning('HRF duration mismatch due to sampling: %f sec',hrf_duration-hrf_time(end));
    end
    
    % Setup HRF
    hrf_types = process_nst_glm_fit('get_hrf_types');
    switch hrf_type
        case hrf_types.CANONICAL  
            hrf = process_nst_glm_fit('cpt_hrf_canonical',hrf_time); 
        case hrf_types.GAMMA 
            hrf = process_nst_glm_fit('cpt_hrf_gamma',hrf_time);
        case hrf_types.BOXCAR 
            hrf = process_nst_glm_fit('cpt_hrf_boxcar',hrf_time);
        otherwise
            bst_error('Unknown hrf_type');
            return;
    end
    
    if size(hrf, 2) ~= 1
        hrf = hrf'; % ensure column vector
    end
    
    %% Make stimulus-induced design matrix    
    X_event = nst_make_event_regressors(events, hrf, time);
    names = {events.label};
    
    
    model.n_roi=n_roi+length(names);
    model.X= [model.X X_event];
    for i_name=1:length(names)
        model.reg_names{end+1}=names{i_name};
    end    

end


function model=nst_glm_add_constant_regressors(model)
    
    model.X= [model.X ones(model.ntime,1)];
    model.reg_names{end+1}='Constant';

end

function model=nst_glm_add_ext_input_regressors(model, sInput_ext,hb_types)
    if ~isempty(sInput_ext) && ~isempty(sInput_ext.FileName)
        
        DataMatExt = in_bst_data(sInput_ext.FileName);
        ChannelExt = in_bst_channel(sInput_ext.ChannelFile);
        dt = diff(DataMat.Time(1:2));
        if length(DataMat.Time) ~= length(DataMatExt.Time) || ...
                ~all(abs(DataMatExt.Time - DataMat.Time) <= dt/100)
            error('Time of external measure is not consistent with data time.');
        end
        
        if length(hb_types)==1
            hb_type=hb_types;
        else
            hb_type = [];
        end    
        
        if ~isempty(hb_type)
            hb_chans = strcmp({ChannelExt.Channel.Group}, hb_type);
            extra_regs = DataMatExt.F(hb_chans, :)';
            extra_reg_names = {ChannelExt.Channel(hb_chans).Name};
        else
            extra_regs = DataMatExt.F';
            extra_reg_names = {ChannelExt.Channel.Name};
        end

         model.X        = [model.X extra_regs];
         model.reg_names= [model.reg_names extra_reg_names];
    end
end


function model=nst_glm_add_DCT_regressors(model,frequences,names)
        
    [cmat,band_indexes] = nst_math_build_basis_dct(model.ntime , model.fs, frequences);
    
    % Add the regressor at the end of the design matrix
    model.X= [model.X cmat]; 
    
    % Add names into the regressor list
    
    for i_name=1:length(names)
        dct_index=band_indexes{i_name};
        for i_index=1:length(dct_index)
           model.reg_names{end+1}= sprintf('%s%i',names{i_name},dct_index(i_index)); 
        end    
    end
    
end
