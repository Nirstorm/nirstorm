function [WData, SData, OPTIONS] = be_discrete_wavelet_transform(Data,OPTIONS,varargin)
% BE_DISCRETE_WAVELET_TRANSFORM computes time-frequency decomposition
%
%   INPUTS:
%       -   Data    :   data matrix (Nsensors x Nsamples)
%       -   OPTIONS :   options structure
%
%   OUTPUTS:
%       -   WData   :   wavelet transform
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: LATIS, 2012
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

        
switch (OPTIONS.wavelet.type)

    case 'RDW'                  
        if OPTIONS.wavelet.vanish_moments == 0
        filtre = MakeONFilter('Haar');
        else
        filtre = MakeONFilter('Daubechies',2*OPTIONS.wavelet.vanish_moments+2);
        end
        [Ns,No] = size(Data);
        Nj = fix(log2(No));
        Noff = min(Nj-1,3);
        WData = zeros(size(Data));
        for i = 1:Ns;
            WData(i,:) = FWT_PO(Data(i,:),Noff,filtre);
        end
        WData(:,1:No/2^(Nj-Noff)) = 0.0;
        SData(:,1:No/2^(Nj-Noff)) = WData(:,1:No/2^(Nj-Noff));
        OPTIONS.automatic.scales(1,:) = 1:Nj-Noff;
        OPTIONS.automatic.scales(2,:) = 3./2.^(1:Nj-Noff)/4;
        
        if OPTIONS.optional.verbose
        disp([OPTIONS.mandatory.pipeline ', Wavelab850 wavelet transform (' OPTIONS.wavelet.type num2str(OPTIONS.wavelet.vanish_moments) '):']);
             for j=1:size(OPTIONS.automatic.scales,2)
                fprintf('\tj = %d,\tfc = %4.3f*fs,\t(%4.1f Hz)\n',OPTIONS.automatic.scales(1,j),OPTIONS.automatic.scales(2,j),OPTIONS.automatic.scales(2,j)*OPTIONS.automatic.sampling_rate);
             end
        end

    case 'rdw'                  
        if OPTIONS.wavelet.vanish_moments > 7
        filter = 'rdw0';
        else
        filter = ['rdw' num2str(OPTIONS.wavelet.vanish_moments)];
        end
        [Ns,No] = size(Data);
        Nj    = fix(log2(No));
        Njs   = max(Nj-3,1);
        SData = zeros(size(Data));
        WData = zeros(size(Data));
        [WData, info ] = be_dwanalysis( Data, Njs, filter );
        WData(:,1:No/2^Njs) = 0.0;
        SData(:,1:No/2^Njs) = WData(:,1:No/2^Njs);
        OPTIONS.automatic.scales(1,:) = 1:Njs;
        OPTIONS.automatic.scales(2,:) = 3./2.^(1:Njs)/4;
        
        if OPTIONS.optional.verbose
        disp([OPTIONS.mandatory.pipeline ', sdw-rdw wavelet transform (' OPTIONS.wavelet.type num2str(OPTIONS.wavelet.vanish_moments) '):']);
             for j=1:size(OPTIONS.automatic.scales,2)
                fprintf('\tj = %d,\tfc = %4.3f*fs,\t(%4.1f Hz)\n',OPTIONS.automatic.scales(1,j),OPTIONS.automatic.scales(2,j),OPTIONS.automatic.scales(2,j)*OPTIONS.automatic.sampling_rate);
             end
        end
       
end       
return