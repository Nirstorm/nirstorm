function results = nst_misc_unpack_glm_result(results, model,method_name,mask)
   

    results.B(:,mask)=model.B;
    if strcmp(method_name,'OLS_prewhitening')
        results.covB(:,:,mask)=model.covB;
    else
        results.covB(:,:,mask)=repmat(model.covB,1,1,length(mask));
    end

    results.residuals(:,mask)=model.residuals;
    results.mse_residuals(:,mask)=model.mse_residuals;
    results.dfe(mask) = model.dfe;

end

