function [dw,OPTIONS] = be_denoise_csoft(w,info_extension,OPTIONS)
% soft schrinkage in the complex case.
%
%   INPUTS: 
%       -   w       :   wavelet transform
%       -   OPTIONS :   options structure
%
%   OUTPUTS: 
%       -   dw      :   denoised transform
%       -   OPTIONS :   options structure
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: JM Lina, 2012, jan 1st
%
%% ==============================================
% License 
%
% BEst is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BEst is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BEst. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------   


% The data size without zero padding added to reach the dyadic length
RealDataSize2Den    = info_extension.end-info_extension.start+1;
% The number of zero padding added to reach the dyadic length
NbZeroPad           = size(w,2)-RealDataSize2Den;

N_level = size(OPTIONS.automatic.scales,2);
[Nc,Nw] = size(w);
dw = w;
killed =[];
    for channel = 1:Nc
        W1 = squeeze(w(channel,Nw/2+1:end));
        
        % Index to take only the real data without zero padding
        i_NoPad_str = round(NbZeroPad/4)+1;
        i_NoPad_end = round(size(W1,2)-NbZeroPad/4);
        W1_NoPad    = W1(:,i_NoPad_str:i_NoPad_end);
        
        sigma_R  = median(abs(real(W1_NoPad)))/0.6745;
        sigma_I  = median(abs(imag(W1_NoPad)))/0.6745;
        thresh   = sqrt(2*log(RealDataSize2Den))*sqrt(sigma_R^2+sigma_I^2);
        
        for iv = 1:N_level
            wi = squeeze(w(channel, 1+Nw/2^iv:Nw/2^(iv-1)));
            d  = (abs(wi) - repmat(thresh,1,length(wi)));
            k = 1./abs(wi);
            k(d<=0) = 0.00;
            killed(channel,iv) = 100*(sum(d<=0)-sum(abs(wi)<eps))/length(wi);
            wi = k.*d.*wi;
            dw(channel,1+Nw/2^iv:2*Nw/2^iv)= wi;
        end
    end
    OPTIONS.automatic.scales = [OPTIONS.automatic.scales ; mean(killed,1)];
    if OPTIONS.optional.verbose
        disp([OPTIONS.mandatory.pipeline ', wavelet denoising:']);
        j=1; fprintf(' j=%d (%d%% to 0),',j,fix(OPTIONS.automatic.scales(3,j)));
        for j=2:size(OPTIONS.automatic.scales,2)
            fprintf(' j=%d (%d%% to 0)',j,fix(OPTIONS.automatic.scales(3,j)));
            if mod(j,3)==0, fprintf('\n'); else fprintf(','); end;
        end
        fprintf('\n');
    end
end