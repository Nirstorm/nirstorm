function fitted_model = nst_glm_fit(model, y, filter_type, filter_param,method,varagin)
%NST_GLM_MODEL_FIT Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4
        method='OLS_precoloring';
    end    

    switch method
        case 'OLS_precoloring'
            [B, covB, dfe, residuals, mse_residuals] = OLS_precoloring_fit(model, y, filter_type, filter_param);
        case 'OLS_prewhitening'
            %if one wants to extend to AR(p) then pass p as varagin{1}
            [B, covB, dfe, residuals, mse_residuals] = OLS_AR1_fit(model, y, filter_type, filter_param);
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


function [B, covB, dfe, residuals, mse_residuals] = OLS_precoloring_fit(model,y, filter_type, filter_param)
    
    % Low pass filter, applied to input data
    lpf = full(process_nst_glm_fit('lpf_hrf',model.hrf, size(y, 1)));
    % Convert to full mat since operation on sparse matrices are
    % single-threaded only but full 
    y_filtered = lpf * y;
    
    if any(~isfinite(y_filtered))
        warning('Found non-finite values in filtered input data.');
    end

    % Band-pass filtering of the design matrix
    model = nst_glm_apply_filter(model,filter_type, filter_param );
    model = nst_glm_apply_filter(model,'lpf', lpf );
    X_filtered=model.X;
    
    %% Solve y = X*B using SVD as in SPM
    [B,proj_X] = nst_glm_fit_B(model,y_filtered, 'SVD');
    
    
    %% For stat afterwards
    res_form_mat = eye(size(lpf)) - X_filtered * proj_X;
    lpf_lpf_T = full(lpf * lpf');
    RV = res_form_mat * lpf_lpf_T;
    trRV = sum(diag(RV));
    RVRVt = RV .* RV';
    trRVRV = sum(RVRVt(:)); % faster than sum(diag(RV * RV));
    dfe = trRV^2 / trRVRV;
    
    fit = X_filtered * B;
    residuals = y_filtered - fit;     
    mse_residuals = var(residuals) * (size(y,1)-1) / trRV;
        
%     figure();plot(y_filtered(:, 161)/max(y_filtered(:,161)), 'b', 'LineWidth', 2); hold on; plot(X_filtered);
%     figure(); hold on; plot(y_filtered(:,161), 'b'); plot(residual(:,161), 'g'); plot(fit(:,161), 'r');
%   
    pXS = proj_X * lpf;
    covB = pXS * pXS';

    if 0 % for test when running GLMTest.test_cortical_simulation
       
        % activ pos hbO: 4708
        % inactiv pos hbo: 4977
        
        poi_inact = 4955;
        
        figure(); hold on;
        plot(y(:,poi_inact), 'k');
        plot(residuals(:,poi_inact), 'g');
        plot(fit(:,poi_inact), 'r');
        mse_residuals_inact = mse_residuals(poi_inact);
        fprintf('MSE_inact = %e\n', mse_residuals_inact);
        t_stat_inact = B(poi_inact) / sqrt(covB(1,1,poi_inact));
        fprintf('tstat_inact = %1.3f\n', t_stat_inact);
        p_val_inact = process_test_parametric2('ComputePvalues', t_stat_inact, dfe, 't', ...
                                               'one+');
        fprintf('p_val_inact = %1.3f\n', p_val_inact);
        
    end
end


function [B_out, covB_out, dfe_out, residuals_out, mse_residuals_out] = OLS_AR1_fit(model,y,filter_type, filter_param)
    
    n_chan=size(y,2);
    n_cond=size(model.X,2);
    n_time=size(y,1);
    
    B_out=zeros(n_cond,n_chan);
    
    covB_out=zeros(n_cond,n_cond,n_chan);
    dfe_out=zeros(1,n_chan);
    residuals_out=zeros(n_time,n_chan);
    mse_residuals_out=zeros(1,n_chan);
    
    % high-pass filtering of the design matrix
     model = nst_glm_apply_filter(model,filter_type, filter_param );
    [B_init,proj_X] = nst_glm_fit_B(model,y, 'SVD');
    bst_progress('start', 'GLM - Pre-whitenning ' , 'Fitting the GLM', 1, n_chan);
    
    for i_chan=1:n_chan
        % Solve B for chan i_chan. We need to solve B for each channel
        % speratly as we are fitting one AR model per channel.
        tic
        fmodel=model; % keep a track of the model so it doesn't change
        
        y=y(:,i_chan);
        SX=fmodel.X;
        SY=y;
       
        B=B_init(:,i_chan);            
        % Estimate the AR(1) processe on the residual
        res=SY-SX*B; 
        W=nst_math_fit_AR(res',1);

        % Compute the filtering matrix 
        a1 = -W(2);
        ka = sparse( ( eye(n_time) - diag( ones(1,n_time-1)*a1,-1) ))^(-1); 
        Va = ka* ka';
        
        % Compute corresponding low-pass filter
        S= full( inv(Va) )^(0.5); % can't use power .5 on sparse matrix
            
        % Apply the filtering to the data and the deisgn matrix
        SY = S*SY;
        fmodel = nst_glm_apply_filter(fmodel,'lpf', S );


        % Compute B for Sy = SXB + Se following an iid normal distribution
        [B,proj_X] = nst_glm_fit_B(fmodel,SY, 'SVD');
        SX=fmodel.X;
        Va=(inv(S))^2;        
        % Compute stat afterwards
       
        pSX=pinv(SX);
        R = eye(n_time) - SX * pSX;
        
        S_Va_S_t = S * Va * S';
        RV = R * S_Va_S_t;
        trRV = sum(diag(RV));
        RVRVt = RV .* RV';
        trRVRV = sum(RVRVt(:)); % faster than sum(diag(RV * RV));
        
        fit = SX * B;
        residuals = SY - fit;
        sigma2=var(residuals);

        dfe = trRV^2 / trRVRV;
        mse_residuals =  sigma2* (size(y,1)-1)/trRV;
        covB = pSX * pSX';

        
        % Save stats 
        B_out(:,i_chan)=B;
        covB_out(:,:,i_chan)=covB;
        dfe_out(i_chan)=dfe;
        residuals_out(:,i_chan)=residuals;
        mse_residuals_out(i_chan)=mse_residuals;
        e = toc;
        disp( [ '#' num2str(i_chan) ' analized in ' num2str(iter) ' iteration (' num2str(e) ' sec)']) 
        bst_progress('inc', 1); 
    end
    
    bst_progress('stop');

end
