function W = nst_math_fit_AR(x,M)
    % Autoregressive all-pole model parameters using Yule-Walker method
    % TODO : do not use signal processing toolbox 
    
    W = aryule(x,M);
    
end