function model = nst_glm_initialize_model(time) 
    model=struct('X',[],'reg_names',{},'n_roi',0,'time',[],'hrf',[],'accept_filter',[], 'ntime',0,'fs',0);
    
    %inizialize the structure
    model(1).time=time;
    model(1).ntime=length(time); 
    model(1).fs=1/(time(2)-time(1));
end

