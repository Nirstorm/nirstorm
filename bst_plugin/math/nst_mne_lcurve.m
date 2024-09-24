function J = nst_mne_lcurve(HM,OPTIONS)
% nst_mne_lcurve - this function solve the inverse probleme using a l-curve
% approach
% Input: HM - struct
%        | - HM.Gain : Gain matrix 
%        OPTIONS - struct
%        | - OPTIONS.Data : data matrix (nChannel x nTimes)
%        | - OPTIONS.DataTime: Corresponding time (1xnTimes)
%        | - OPTIONS.BaselineSegment: Segment used for baseline [start,end]
%        | - OPTIONS.NoiseCov_recompute: true if using the baseline to
%        estimate the noise covariance, if false, the noise covariance is
%        identity
%        | - OPTIONS.depth_weigth_MNE: depth-weighting factor between 0 and 1
% Output: J - reconstructed time-course on the cortex (nVertex x nTimes)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Edouard Delaire, 2024

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
        Sigma_s = diag(power(diag(G'*inv(Sigma_d)*G),-p)); 
    else
        Sigma_s = diag(power(diag(G'*G),-p)); 
    end
    

    % Pre-compute matrix
    GSG = G * Sigma_s * G';
    SG  = Sigma_s * G';

    % Note W*S*G: W = Sigma_s^-0.5 so W*Sigma_s = Sigma_s^0.5
    wSG = sqrt(Sigma_s) * G';

    % Parameter for l-curve
    param1  = [0.1:0.1:1 1:5:100 100:100:1000]; 

    % Scale alpha using trace(G*G')./trace(W'*W)  
    scale   = trace(G*G')./ trace(inv(Sigma_s));       
    alpha   = param1.*scale;

    % Pre-compute data decomposition
    [U,S]   = svd(M,'econ'); 


    Fit     = zeros(1,length(alpha));
    Prior   = zeros(1,length(alpha));

    bst_progress('start', 'wMNE, solving MNE by L-curve ... ' , 'Solving MNE by L-curve ... ', 1, length(param1));
    for iAlpha = 1:length(alpha)
        
        inv_matrix = inv( GSG  + alpha(iAlpha) * Sigma_d );
        
        residual_kernal = eye(size(M,1)) - GSG * inv_matrix;
        wKernel = wSG*inv_matrix;


        % Estimate the corresponding norm
        R = qr(residual_kernal*U);
        Fit(iAlpha)     = norm(R*S);

        R = qr(wKernel*U);
        Prior(iAlpha)   = norm(R*S);

    
        bst_progress('inc', 1); 
    end

    % Fid alpha optimal based on l-curve
    [~,Index] = min(Fit/max(Fit)+Prior/max(Prior)); 
    
    Kermel = SG * inv( GSG  + alpha(Index) * Sigma_d );
    J = Kermel*M; 


    bst_progress('text', 'wMNE, solving MNE by L-curve ... done');
end
