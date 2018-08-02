function [B,covB,dfe]=nst_ols_fit(y,X)
    B= pinv(X'*X)* X'*y; % or B=X\y; 
    
    residual=y - X*B ;     
    covE=cov(residual);
    
    for i=1:size(residual,2)     
        covB(:,:,i)=covE(i,i) * pinv(transpose(X)*X);
    end
    dfe = size(y,1) - rank(X);
    
end