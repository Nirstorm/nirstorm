function J = nst_mne_lcurve_MAP(HM,OPTIONS)
% nst_mne_lcurve - this function solve the inverse probleme using a l-curve
% approach in the MAP formalism. This approach is inefficient and should
% not be used. Consider using nst_mne_lcurve instead.
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
        Pd    =   inv(diag(diag(real(cov(baseline')))));
    else
        Pd    =   eye(size(M,1));  
    end

    p = OPTIONS.depth_weigth_MNE;
    if OPTIONS.NoiseCov_recompute
        Ps = diag(power(diag(G'*Pd*G),p)); 
    else
        Ps = diag(power(diag(G'*G),p)); 
    end
    W = sqrt(Ps);

    % Parameter for l-curve
    param1  = [0.1:0.1:1 1:5:100 100:100:1000]; %the param1 list we tested in wMNE_org

    % Scale alpha using trace(G*G')./trace(W'*W)  
    scale   = trace(G*G')./ trace(Ps) ;       
    alpha   = param1.*scale;
    
    % Pre-compute
    GPG = G'*Pd*G;
    GP  = G'*Pd;
    Fit     = zeros(1,length(param1));
    Prior   = zeros(1,length(param1));

    bst_progress('start', 'wMNE, solving MNE by L-curve ... ' , 'Solving MNE-MAP by L-curve ... ', 1, length(param1));
    for iAlpha = 1:length(alpha)
        Kernel = inv( GPG  +  alpha(iAlpha) * Ps)*GP;
        J =  Kernel * M; 

        Fit(iAlpha)     = norm(M-G*J);  % Define Fit as a function of alpha
        Prior(iAlpha)   = norm(W*J);      % Define Prior as a function of alpha
    
        bst_progress('inc', 1); 
    end

    % Fid alpha optimal based on l-curve
    [~,Index] = min(Fit/max(Fit)+Prior/max(Prior)); 

    Kernel = inv( GPG  +  alpha(Index) * Ps)*GP;
    J =  Kernel * M; 

    bst_progress('text', 'wMNE, solving MNE by L-curve ... done');

end