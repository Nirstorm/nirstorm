function DEF = be_wmem_pipelineoptions()

        % clustering
        DEF.clustering.clusters_type        = 'wfdr';
        DEF.clustering.MSP_scores_threshold = 'fdr';
%        DEF.clustering.neighborhood_order   = 4;
        
        % model
        DEF.model.alpha_threshold       = 0.10;
        DEF.model.active_mean_method    = 2;
        DEF.model.alpha_method          = 3;
        
        % wavelet processing
        DEF.wavelet.type                = 'rdw';
        DEF.wavelet.vanish_moments      = 4;
        DEF.wavelet.shrinkage           = 1;
        DEF.wavelet.selected_scales     = 0;
        DEF.wavelet.verbose             = 0;
        DEF.wavelet.single_box          = 0;
        
        % automatic
        DEF.automatic.selected_samples  = [];
        DEF.automatic.selected_jk       = [];
        DEF.automatic.selected_values   = [];
        DEF.automatic.Mod_in_boxes      = [];
        DEF.automatic.scales            = [];
        
        % solver
        DEF.solver.spatial_smoothing    = 0.6;
        DEF.solver.Optim_method         = 'fminunc';
        DEF.solver.NoiseCov_method      = 5;
        
        % optional
        DEF.optional.normalization      = 'fixed'; 

return