function model = nst_glm_apply_filter(model,filter_name, varargin )
%NST_GLM_APPLY_FILTER Apply filter on the regressor of the column matrix
%but only apply when it's ok to apply it (for exemple, we don't apply a
%high-pass filter on the constant regressor). Condition about what filter
%can be applied can be found in model.accept_filter. For a channel i, if
% model.accept_filter = 0 then no filter can be applied
% model.accept_filter = 1 then high-pass can be applied
% model.accept_filter = 2 then low-pass can be applied
% model.accept_filter = 3 then  both can be applied


    switch filter_name
        
        case 'lpf'
            % Applied the GLM low pass filter;
            % Usage : 
            % lpf = full(lpf_hrf(hrf, size(y, 1)));
            % model = nst_glm_apply_filter(model,'lpf', lpf );

            lpf=varargin{1};
            
            ind=find(model.accept_filter >= 2);

            model.X(:,ind)= lpf*model.X(:,ind);
            model.X(:,ind) = model.X(:,ind) ./ max(abs(model.X(:,ind)));

        case 'IIR_highpass'
            % Applied the GLM IIR high pass filter;
            % Usage  model = nst_glm_apply_filter(model,'IIR_highpass', hpf_low_cutoff );
            cutoff=varargin{1};
            if cutoff > 0 % Security check :)
                ind=find(model.accept_filter == 1 | model.accept_filter == 3);
                model.X(:,ind)= process_nst_iir_filter('Compute', model.X(:,ind), model.fs, ...
                                                       'highpass', cutoff, 0, 2, 0);
            end
        case 'IIR_bp'
            % Applied the GLM IIR high pass filter;
            % Usage  model = nst_glm_apply_filter(model,'IIR_highpass', hpf_low_cutoff );
            low_cutoff = varargin{1};
            high_cutoff = varargin{2};
            if nargin < 5
                order = 3;
            else
               order = 3; varargin{3};
            end

            ind=find(model.accept_filter == 1 | model.accept_filter == 3);
            model.X(:,ind)= process_nst_iir_filter('Compute', model.X(:,ind), model.fs, ...
                                                       'bandpass', low_cutoff, high_cutoff, order, 0);
        case 'FIR_bp'
            low_cutoff = varargin{1};
            high_cutoff = varargin{2};
            
            if nargin < 5
                is_relax = 1;
            else
                is_relax = varargin{3};
            end
            if nargin < 6
                TranBand = min(0.5*low_cutoff, 0.2*high_cutoff); % completely arbritary
            else
                TranBand = varargin{4};
            end

            ind=find(model.accept_filter == 1 | model.accept_filter == 3);
            
            [x, FiltSpec, Messages] = process_bandpass('Compute', model.X(:,ind)',  model.fs, low_cutoff, high_cutoff, [], 0,is_relax , TranBand);
            model.X(:,ind)= x;
   
        case 'DCT_filter'
            % Applied the GLM IIR DCT filter;
            % Usage  model = nst_glm_apply_filter(model,'DCT', low_period_cuttoff );
            cutoff=varargin{1};
            
            model_dct= nst_glm_initialize_model(model.time);
            model_dct=nst_glm_add_regressors(model_dct,'constant');
            model_dct=nst_glm_add_regressors(model_dct,'linear');  

            if cutoff > 0
                model_dct=nst_glm_add_regressors(model_dct,'DCT',[1/model.time(end) 1/cutoff],{'LFO'});  
            end
            
            ind=find(model.accept_filter == 1 | model.accept_filter == 3);
            [B,proj_X] = nst_glm_fit_B(model_dct,model.X(:,ind), 'SVD');
            
            model.X(:,ind) = model.X(:,ind) - model_dct.X*B;
    end    

end

