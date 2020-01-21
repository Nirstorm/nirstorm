function [out]= be_sdwsynthesis(in_SDW, Nb_Level, filter)
% Complex Daubechies wavelet inverse transform
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
%  Authors: JM Lina, Benoit Decarie, 2012, jan 1st
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


    filter_list = {'sdw0','sdw2','sdw6','sdw8'};
    if ~ismember(filter,filter_list)
        disp('!! invalid complex filter: we use sdw2')
        filter = 'sdw2';
    end
    [analF analG H0 G0] = be_makeqfb(filter);

    % data are 1xN or NcxN with N = power of 2
    flip = 0;
    if size(in_SDW,2)==1
        in_SDW = in_SDW'; flip = 1;
    end
    
    
    [Nb_line extend trial] = size(in_SDW);  % Dimensions of the input.
    dim_H0 = size(H0,2);            % Get number of columns.
    Jcase    = (dim_H0-2)/2;        % Number of zero momentum of the HP filter.

    % Division des parties reelles et imaginaires.
    in_re = real(in_SDW); in_im = imag(in_SDW);
    H_re = real(H0);    H_im = imag(H0);
    G_re = real(G0);    G_im = imag(G0);
    
    scale = bitshift(1,Nb_Level);   % (1 << Nb_Level); 
    n = extend/scale;                    % taille de V
    out_re = in_re; out_im = in_im;         % Initialisation.
    j = n;
    for ii = 1:Nb_Level
        n = bitshift(n,1); 
        
        temp_a  = fix_vector(out_im, n);    
        temp_b  = fix_vector(in_im(:,j+1:n,:), n);
        temp_im = be_convsynthese(temp_a, temp_b, H_re, G_re, n, Jcase);
        temp_re = be_convsynthese(temp_a, temp_b, H_im, G_im, n, Jcase);

        temp_a  = fix_vector(out_re, n);    
        temp_b  = fix_vector(in_re(:,j+1:n,:), n);
        out_im  = be_convsynthese(temp_a, temp_b, H_im, G_im, n, Jcase);
        out_re  = be_convsynthese(temp_a, temp_b, H_re, G_re, n, Jcase);

        out_re(:,1:n,:) = out_re - temp_re; 
        out_im(:,1:n,:) = out_im + temp_im;

        j = bitshift(j,1);  % left shifting "j <<= 1".
    end
    
    % On regroupe les parties reelles et imaginaires.
    out = out_re + 1i*out_im;
    if flip , out = out'; end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = fix_vector(in, num)
%START
    sz = size(in);
    if (size(sz,2) == 2)
        sz(3) = 1;
    end
    
    out = zeros(sz(1),num,sz(3)); % Preallocation.
    j = 1; % Initialisation.
    for i =1:2:num
        out(:,i,:) = in(:,j,:);
        j = j+1;
    end
end