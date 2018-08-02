function [B,covB,dfe]=nst_ar_irls_fit(y,X,pmax)
	stat=nirs.math.ar_irls(y,X, pmax );
    
    B=stat.beta;
    
    covB=zeros( size(stat.covb,1), size(stat.covb,2),size(stat.covb,3));
    for i=1:size(stat.covb,3)
       covB(:,:,i)= stat.covb(:,:,i,i);
    end    
    dfe=stat.dfe;
end