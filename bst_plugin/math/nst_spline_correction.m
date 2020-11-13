function nirs_sig_corr = nst_spline_correction(nirs_sig, t, mvtWin_samples)
%% Motion correction algorithm based on spline interpolation 
%  Scholkmann, F., Spichtig, S., Muehlemann, T., & Wolf, M. (2010). 
%  How to detect and reduce movement artifacts in near-infrared imaging 
%  using moving standard deviation and spline interpolation. 
%  Physiological measurement, 31(5), 649–662. 
%  https://doi.org/10.1088/0967-3334/31/5/004
%
% Args
%    - nirs_sig: matrix of double, size: time x nb_channels
%        nirs signals to be corrected
%    - t: array of doublee, size: time x 1
%        time axis of the given nirs signal
%    - mvtWin_samples: matrix of double, size: nb_mvt_events x 2
%        windows of marked movements (starting and ending sample index)
%
% Output:
%    - nirs_sig_corr: matrix of double, size: time x nb_channels
%        Motion-corrected nirs signals
%
param.spline.smooth=0.99;

nSamples = size(nirs_sig, 1);
nChan = size(nirs_sig, 2);
nMvt = size(mvtWin_samples, 1);
d = nirs_sig;
%====================================================================
% define normal and mvt segments
%======================================================================


for iM=1:nMvt
    if iM==1
        normalWin_samples(iM,:)=[1 mvtWin_samples(iM,1)-1];
    else
        normalWin_samples(iM,:)=[mvtWin_samples(iM-1,2)+1 mvtWin_samples(iM,1)-1];
    end
    
    if iM==nMvt
        normalWin_samples(iM+1,:)=[mvtWin_samples(iM,2)+1 nSamples];
    end
end
clear iM
nNormalSeg=size(normalWin_samples,1);

%==========================================================================
% Preprocessing MC correction: all channels
%==========================================================================
% F Scholkmann; Physiol.Meas.31 2010

% =============================================================
% segment data
% =============================================================
for iM=1:nMvt
    MA_seg{iM}= d(mvtWin_samples(iM,1):mvtWin_samples(iM,2),:);
end

for iN=1:nNormalSeg
    normal_seg{iN}= d(normalWin_samples(iN,1):normalWin_samples(iN,2),:);
    normal_seg_mean(iN,:)=mean(normal_seg{iN},1);
end
normal_seg_mean_global=mean(normal_seg_mean);

% =============================================================
% estimate the mvt artefact: spline with smooth parameter
% =============================================================
for iM=1:nMvt
    for iChan=1:nChan
        interp_seg{iM}(:,iChan) = fnval(...
            csaps(...
            t(mvtWin_samples(iM,1):mvtWin_samples(iM,2)),...
            MA_seg{iM}(:,iChan),...
            param.spline.smooth),...
            t(mvtWin_samples(iM,1):mvtWin_samples(iM,2)));
    end
end
clear iChan


% =============================================================
% substract the estimated spline from the MA periods
% =============================================================
for iM=1:nMvt
    denoised_seg{iM} = MA_seg{iM} - interp_seg{iM};
    % off=repmat(mean(d([mvtWin_samples(iM,1)-1 mvtWin_samples(iM,2)+1],:)),size(MA_seg{iM},1),1);
end

% =============================================================
% shift all normal segments and replace them in d:
% the method here differ a little from the one proposed by  F Scholkmann
% =============================================================

% -------------------------------------------------------------------------
% % Solution 1: fit to the last sample of the previous segment
% -------------------------------------------------------------------------
%             for iN=2:nNormalSeg
%                 nS_iSeg=size(normal_seg{iN},1);
%                 d(normalWin_samples(iN,1):normalWin_samples(iN,2),:)=normal_seg{iN}-...
%                                                                      repmat(normal_seg_mean(iN,:),nS_iSeg,1)+...
%                                                                      repmat(d(normalWin_samples(iN-1,2),:),nS_iSeg,1);
%             end

% -------------------------------------------------------------------------
% % Solution 2: weighted shift of segments according to the
% number of samples before and after the transition
% -------------------------------------------------------------------------
for iN=1:nNormalSeg-1
    
    nS_iSeg_1=size([1:normalWin_samples(iN,2)]',1);
    nS_iSeg_2=size(d(normalWin_samples(iN,2)+1:end,:),1);
    
    nS_iSeg_1_2=nS_iSeg_1+nS_iSeg_2;
    weight_iSeg_1=(nS_iSeg_1/nS_iSeg_1_2)/0.5;
    weight_iSeg_2=(nS_iSeg_2/nS_iSeg_1_2)/0.5;
    
    val_seg_1=d(normalWin_samples(iN,2),:);
    val_seg_2=d(normalWin_samples(iN+1,1),:);
    shift=( val_seg_2-val_seg_1)/2;
    
    d(1:normalWin_samples(iN,2),:)=d(1:normalWin_samples(iN,2),:)+...
        repmat(weight_iSeg_2*shift,normalWin_samples(iN,2),1);
    
    d(normalWin_samples(iN,2)+1:end,:)=  d(normalWin_samples(iN,2)+1:end,:) -...
        repmat(weight_iSeg_1*shift,size(d(normalWin_samples(iN,2)+1:end,:),1),1);
end



% =============================================================
% add a  constant offset  into the MA residuals to fit previous
% and next normal segments ( interp function is used so if we
% choose another method to shift segments we can add linear
% offset instead)
% =============================================================
for iM=1:nMvt
    off=[];
    for iC=1:nChan
        off(:,iC)=interp1([mvtWin_samples(iM,1)-1;mvtWin_samples(iM,2)+1],...
            [d(mvtWin_samples(iM,1)-1,iC);d(mvtWin_samples(iM,2)+1,iC)],...
            [mvtWin_samples(iM,1)-1:mvtWin_samples(iM,2)+1]',...
            'linear');
    end
    denoised_seg{iM} = denoised_seg{iM} + off(2:end-1,:);
end

% =============================================================
% replace MA segments in d
% =============================================================
for iM=1:nMvt
    d(mvtWin_samples(iM,1):mvtWin_samples(iM,2),:)=denoised_seg{iM};
end

% =============================================================
% offset <0 channels
% =============================================================
% [d ~]=mfip_correct_negChan(d,ml);

nirs_sig_corr = d;

end

