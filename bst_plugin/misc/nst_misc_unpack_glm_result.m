function [B,covB,dfe,residuals,mse_residuals] = nst_misc_unpack_glm_result(model,method_name,surface_data,nb_regressors,n_voxel,mask)
    if surface_data
        B=zeros(nb_regressors,n_voxel);
        B(:,mask)=model.B;
        
        residuals=zeros( size(model.residuals,1),n_voxel);
        residuals(:,mask)=model.residuals;

        mse_residuals=zeros(1,n_voxel);
        mse_residuals(:,mask)=model.mse_residuals;
        
        dfe=zeros(1,n_voxel);
        dfe(:,mask)=model.dfe;
        
        covB=model.covB;
    else 
        B=model.B;
        residuals=model.residuals;
        mse_residuals=model.mse_residuals;
        
        dfe=model.dfe;
        covB=model.covB;     
    end
    
    if strcmp(method_name,'OLS_prewhitening') && surface_data
        covB=zeros(nb_regressors,nb_regressors,n_voxel);
        covB(:,:,mask)=model.covB;

        dfe=zeros(1,n_voxel);
        dfe(mask)=model.dfe;        
    end
    
    if strcmp(method_name,'OLS_precoloring') && surface_data
        covB=model.covB; % This is weird
    end    
end

