function [out_RDW, info] = be_rdwanalysis(in_C, Nb_Levels, filter)
% Real Daubechies wavelet transform (using Wavelab850)
%
%   INPUTS: 
%       -   in_C    :   (Nc x Nt data)
%       -   Nb_Level:   (must be < log2(Nt))
%       _   filter  :   (integer)
%
%   OUTPUTS: 
%       -   out_SDW :   (Nc x Nt matrix of wavelet coeff.)
%       -   info    :   (structure that keeps the information of the transform)
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


    % WE ASSUME THE DATA IN 2 DIMENSIONS
    % data are 1xN or NcxN with N = power of 2
    [Nb_line Nb_samples] = size(in_C);
    flip = 0;
    if Nb_samples == 1
    in_C = in_C'; 
    [Nb_line Nb_samples] = size(in_C);
    flip = 1;
    end
    % dyadic extension of rhe data (if necessary)
    [in_C, info] = be_extended_dyadic( in_C );
    % we construct the filter
    nfilter = MakeONFilter('Daubechies',filter);
    % we verify the number of levels
    N = ceil(log2(Nb_samples));
    Nb_Levels = max(0,min(N-3,Nb_Levels));
    % initialisation
    out_RDW = zeros(size(in_C));
    % loop over the channels:
    for chan = 1:Nb_line
           sig = in_C(chan,:);
           out_RDW(chan,:) = FWT_PO(sig,N-Nb_Levels,nfilter);
    end
    
    if flip , out_RDW = out_RDW'; end
    % The info:
    info.level = Nb_Levels;
    info.filter = filter;
end