function [analF analG synF synG Jcase] = makeQFB(filtre)
% Complex SDW Filter bank construction
% Ref: JM Lina, Journ. Math. of Vision (1995)
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

    load(['be_' filtre]);
    filter_length = 2*filtre.J+2;
    analF = zeros(1,filter_length);
    analG = zeros(1,filter_length);
    synF  = zeros(1,filter_length);
    synG  = zeros(1,filter_length);
    filtre.H
    sg = 1.0;
    Jcase=(length(filtre.H)-2)/2;
    for k = 1:filter_length
       re = real(filtre.H(k));
       im = imag(filtre.H(k));
       analF(k) = filtre.H(k);
       analG(k) =  sg*re-sg*im*1i;
       synF(k)  =  re-im*1i;
       synG(k)  = -sg*re-sg*im*1i;
       sg       = -sg;
    end