function J = nst_mne_lcurve(HM,OPTIONS)

    % selection of the data:
    sample_baseline = be_closest(OPTIONS.BaselineSegment([1 end]), OPTIONS.DataTime);
    baseline        = OPTIONS.Data(:,sample_baseline(1):sample_baseline(2));
    
    sample_data     = be_closest(OPTIONS.TimeSegment([1 end]), OPTIONS.DataTime);
    M               = OPTIONS.Data(:,sample_data(1):sample_data(2));

    % Normalization by using the mean standard deviation of baseline
    
    SD = std(baseline,[],2);
    MSD = mean(SD);

    %Normalize datas, baseline, gain matrix  with the mean std dev
    M           = M./MSD;
    baseline    = baseline./MSD;
    G           = HM.Gain./MSD;
    
    % Compute inverse of noise covariance matrices
    if OPTIONS.NoiseCov_recompute 
        Sigma_d    =   diag(diag(real(cov(baseline'))));
    else
        Sigma_d    =   eye(size(M,1));  
    end

    p = OPTIONS.depth_weigth_MNE;
    if OPTIONS.NoiseCov_recompute
        Ps = diag(power(diag(G'*inv(Sigma_d)*G),p)); 
    else
        Ps = diag(power(diag(G'*G),p)); 
    end
    W = sqrt(Ps);
    Sigma_s = inv(Ps);

    % Pre-compute matrix
    GSG = G * Sigma_s * G';
    SG  = Sigma_s * G';

    
    % Parameter for l-curve
    param1  = [0.1:0.1:1 1:5:100 100:100:1000]; 

    % Scale alpha using trace(G*G')./trace(W'*W)  
    scale   = trace(G*G')./ trace(Ps) ;       
    alpha   = param1.*scale;


    Fit     = zeros(1,length(alpha));
    Prior   = zeros(1,length(alpha));

    bst_progress('start', 'wMNE, solving MNE by L-curve ... ' , 'Solving MNE by L-curve ... ', 1, length(param1));
    for iAlpha = 1:length(alpha)
        
        Kermel = SG * inv( GSG  + alpha(iAlpha) * Sigma_d );
        J = Kermel*M; 

        Fit(iAlpha)     = norm(M-G*J);  % Define Fit as a function of alpha
        Prior(iAlpha)   = norm(W*J);      % Define Prior as a function of alpha
    
        bst_progress('inc', 1); 
    end

    % Fid alpha optimal based on l-curve
    [~,Index] = min(Fit/max(Fit)+Prior/max(Prior)); 
    
    Kermel = SG * inv( GSG  + alpha(Index) * Sigma_d );
    J = Kermel*M; 


    bst_progress('text', 'wMNE, solving MNE by L-curve ... done');
end
