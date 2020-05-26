function fitted_model = nst_glm_fit(model, y, hpf_low_cutoff,method,varagin)
%NST_GLM_MODEL_FIT Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4
        method='OLS_precoloring';
    end    
    dt=model.time(2)-model.time(1);

    switch method
        case 'OLS_precoloring'
            [B, covB, dfe, residuals, mse_residuals] = process_nst_glm_fit('ols_fit',y, dt, model.X, model.hrf, hpf_low_cutoff);
        case 'OLS_prewhitening'
            [B, covB, dfe, residuals, mse_residuals] = process_nst_glm_fit('AR1_ols_fit',Y, dt, model.X, hpf_low_cutoff);
        otherwise
            bst_error('This method is not implemented');
            return
    end
    
    fitted_model=struct('B',[],'covB',[],'dfe',[],'residuals',[],'mse_residuals',[]);
    fitted_model.B=B;
    fitted_model.covB=covB;
    fitted_model.dfe=dfe;
    fitted_model.residuals=residuals;
    fitted_model.mse_residuals=mse_residuals;
end

