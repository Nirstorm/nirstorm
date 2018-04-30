function DEF = be_rmem_pipelineoptions()

        % clustering
        DEF.clustering.clusters_type    = 'blockwise';
        DEF.clustering.MSP_window       = 10;
        DEF.clustering.MSP_scores_threshold = 0;
                
        % wavelet processing
        DEF.wavelet.type                = 'CWT';
        DEF.wavelet.vanish_moments      = 4;
        DEF.wavelet.order               = 10;
        DEF.wavelet.nb_levels           = 128;
        DEF.wavelet.verbose             = 0;
        
        % ridge processing
        DEF.ridges.scalo_threshold      = .95;      
        DEF.ridges.energy_threshold     = .95; 
        DEF.ridges.strength_threshold   = [NaN];
        DEF.ridges.skim_map             = 1;        
        DEF.ridges.frequency_range      = [];
        DEF.ridges.min_duration         = []; %(ms)
        DEF.ridges.cycles_in_window     = 2;
        DEF.ridges.fdrMOD               = {};
        DEF.ridges.method               = 'inhouse';%'LillyOlhede2010';
        
        % automatic
        DEF.automatic.selected_samples  = [];
        DEF.automatic.selected_jk       = [];    
        DEF.automatic.rMEMfiles         = {};
        DEF.automatic.Mod_in_ridges     = {};
        
        % optional
        DEF.optional.normalization      = 'adaptive'; 

        % solver
        DEF.solver.NoiseCov_method      = 2;
        
        % model
        DEF.model.alpha_threshold       = 0.10;
        DEF.model.active_mean_method    = 2;
        DEF.model.alpha_method          = 3;
        
return