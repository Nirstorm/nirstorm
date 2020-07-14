function [B,proj_X] = nst_glm_fit_B(model,y, method)
% nst_glm_model_fit_B - fit the model to compute B using either SVD or the
% pseudo-inverse
    if nargin < 3
        method='SVD';
    end    

    X=model.X;
    switch method
        case 'SVD'
            [X_u, X_s, X_v] = svd(X,0);
            X_diag_s = diag(X_s);
            X_rank_tol =  max(size(X)) * max(abs(X_diag_s)) * eps;
            X_rank =  sum(X_diag_s > X_rank_tol);

            % Compute projector based on X to solve Y = X*B
            proj_X = X_v(:,1:X_rank) * diag(1./X_diag_s(1:X_rank)) * X_u(:,1:X_rank)';
    
            % Fit to data
            B = proj_X * y; 
        case 'pinv'
            proj_X= (X'*X)\X'; % compute inv(X'X)*X'
            B=proj_X*y;
    end          
end

