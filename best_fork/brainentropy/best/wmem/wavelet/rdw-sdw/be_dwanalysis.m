function [ out_SDW, info ] = be_dwanalysis( in_C, Nb_Level, filter )
% Real Daubechies wavelet transform
% 
%   INPUTS: 
%       -   in_C    :   (Nc x Nt data)
%       -   Nb_Level:   (must be < log2(Nt))
%       -   filter  :   ('sdw0, 2, 4, 6 and 8')
%
%   OUTPUTS: 
%       -   out_SDW :   (Nc x Nt matrix of wavelet coeff.)
%       -   info    :   (structure that keeps the information of the transform)
%
% Reference: JM. Lina, Signal Processing with Complex Daubechies Wavelets,
%            Journ. Math. Imaging, 1995
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: JM Lina, Benoit Decarie, Xavier Leturc 2014, july 1st
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

% WE ASSUME THE DATA IN 2 DIMENSIONS
    % data are 1xN or NcxN with N = power of 2
    reel_fil=0;
    [Nb_line Nb_samples] = size(in_C);
    flip = 0;
    if Nb_samples == 1
    in_C = in_C'; 
    [Nb_line Nb_samples] = size(in_C);
    flip = 1;
    end
    % we validate the filter
    be_what_filter_list;
    if ~ismember(filter,filter_list)
        if filter(1)=='r'
            disp('!! invalid real filter: we use rdw2')
            filter = 'rdw2';
        elseif filter(1)=='s'
            disp('!! invalid complex filter: we use sdw2')
            filter = 'sdw2';
        else
            disp('!! invalid filter: we use rdw2')
            filter = 'rdw2';
        end
    end
    if filter(1)=='r' % Cas r�el
        reel_fil=1;
    end
    if reel_fil==1
        
        [H0 G0 synF synG Jcase] = be_makeqfbreal(filter);
    else
        [H0 G0 synF synG Jcase] = be_makeqfb(filter);
    end
  
    % let us compute the SDW transform:
    
%     dim_H0 = size(H0,2);      % Get number of columns.
%     Jcase  = (dim_H0-2)/2;    % Number of zero momentum of the HP filter.
%     Jcase

    % real part and imaginary part are handled separately:
    
    sig=in_C;
    out_re = real(in_C); out_im = imag(in_C);
    H_re = real(H0);    H_im = imag(H0);
    G_re = real(G0);    G_im = imag(G0);
    
    % Preallocation
    %trial (last dimension) = 1 here.
    
    temp_re = zeros(Nb_line, Nb_samples, 1);    
    temp_im = zeros(Nb_line, Nb_samples, 1);
    re = zeros(Nb_line, Nb_samples, 1);
    im = zeros(Nb_line, Nb_samples, 1);

    temp=zeros(Nb_line, Nb_samples, 1);  
    matr=zeros(Nb_line, Nb_samples, 1);
    temp_2=sig;
    m = Nb_samples; % Initialisation.
    n = 0;
    
    if reel_fil == 1
         for i = 1:Nb_Level
            if mod(m,2)==0 
                [temp(:,1:m/2,:) temp(:,m/2+1:m,:)] = ...
                be_convanalysisreal( temp_2, H0, G0, m, Jcase); 
                temp_2=temp(:,1:m/2,:);
                m = bitshift(m,-1);
                n = n+1;
            end
         end
    else
    for i = 1:Nb_Level
        if mod(m,2)==0
            [temp_im(:,1:m/2,:) temp_im(:,m/2+1:m,:)] = ...
                be_convanalysis( out_re, H_im, G_im, m, Jcase);
            [temp_re(:,1:m/2,:) temp_re(:,m/2+1:m,:)] = ...
                be_convanalysis( out_im, H_im, G_im, m, Jcase);
            [re(:,1:m/2,:) re(:,m/2+1:m,:)] = ...
                be_convanalysis( out_re, H_re, G_re, m, Jcase); 
            [im(:,1:m/2,:) im(:,m/2+1:m,:)] = ...
                be_convanalysis( out_im, H_re, G_re, m, Jcase);

            out_re(:,1:m,:) = re(:,1:m,:) - temp_re(:,1:m,:);
            out_im(:,1:m,:) = im(:,1:m,:) + temp_im(:,1:m,:);
       
            [temp(:,1:m/2,:) temp(:,m/2+1:m,:)] = ...
            be_convanalysis( temp_2, H0, G0, m, Jcase); 
            temp_2=temp(:,1:m/2,:);
            m = bitshift(m,-1);
            n = n+1;
        end
    end
    end
    % The complex output:
    if reel_fil==1
        out_SDW = temp;
    else
        out_SDW = out_re + out_im.*(0+1i);
    end

    if flip , out_SDW = out_SDW'; end
    % The info:
    info.level = n;
    info.filter = filter;

