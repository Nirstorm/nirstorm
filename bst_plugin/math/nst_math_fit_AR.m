function W = nst_math_fit_AR(x,M)
    % Autoregressive all-pole model parameters using Yule-Walker equation   
    % W = aryule(x,M);
    
    r = xcorr(x,'biased')'; 
    r(1:length(x)-1) = [];
    
    R=toeplitz(r(1:M));
    
    W=R\(-r(2:M+1));
    W=[1 W'];
end