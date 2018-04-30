function [disl, WData, OPTIONS] = be_get_scalogram(WData, OPTIONS)
%BE_GET_SCALOGRAM computes a normalized scalogram from a time-frequency plane
%
%   INPUTS:
%       -   WData   :   time-frequency plane
%       -   OPTIONS :   options structure
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

% Gets T-F plane dimensions
[nbF, nbS]  = size(WData);
slice       = zeros(nbF, nbS);
lim1        = max(1,be_closest(OPTIONS.ridges.frequency_range(1), OPTIONS.wavelet.freqs_analyzed)-1);
lim2        = min(nbF,be_closest(OPTIONS.ridges.frequency_range(2), OPTIONS.wavelet.freqs_analyzed) + 1);

% Calcul des scalogrammes normalises
sc_mat                      =   OPTIONS.wavelet.scale_analyzed' *  ones(1,nbS);
O                           =   OPTIONS;
O.optional.TimeSegment      =   O.mandatory.DataTime;
O.ridges.strength_threshold =   .5;
O.ridges.skim_map           =   0;

% Normalize the scalogram with the scale
slice(lim1:lim2,:)          = squeeze( WData(lim1:lim2,:) ) ;
slice                       = abs(slice.^2) ./ sc_mat;
clear sc_mat 

% Get local maxima 
[disl]                      =   be_localmaxima( slice, O, [lim1,lim2] );

return



function [ISVALID] = validate_rdg_point(plane, point, crit)

    %%% 1st degree neighborhood %%%
    %
    %   - - - - -
    %   - n n n -
    %   - n p n -
    %   - n n n -
    %   - - - - -
    %
    [N,M] = size(plane);
    y = mod(point,N);
    x = ceil(point/N);
    xs = x : x+1; ys = y-1 : y+1;
    xs(xs>M) = M; ys(ys>N) = N; ys(ys<1)=1;
    COLONNE	= [ys ys];
    LIGNE	= [xs; xs; xs];
    DEG1	= ( LIGNE(:)' - 1)*N + COLONNE;
    %
    %
    %%% END 1st degree neighborhood %%%
    
    
    SCORE = sum( plane( myunique(DEG1,1) ) ) - 1;
    if SCORE >= crit
        ISVALID = true;
        return
    elseif SCORE == 0
        ISVALID = false;
        return
    end
    
    
    %%% 2nd degree neighborhood %%%
    %
    %   n n n n n
    %   n - - - n  
    %   n - p - n
    %   n - - - n
    %   n n n n n
    %
    xs = x : x+2; ys = y-2 : y+2;
    ys(ys<1) = 1; xs(xs>M) = M; ys(ys>N) = N;
    COLONNE	= [ys ys ys ];
    LIGNE	= [xs; xs; xs; xs; xs];
    DEG2	= ( LIGNE(:)' - 1 )*N + COLONNE;
    TRACK = zeros(1,5);
    TRACK(1) = sum( plane( myunique(DEG2([6 11 16]),0) ) ) * plane( DEG1(1) ) * ~( DEG1(1)==DEG2(11) );
	TRACK(2) = sum( plane( myunique(DEG2([10 15 20]),0) ) ) * plane( DEG1(3) ) * ~( DEG1(3)==DEG2(15) );
	TRACK(3) = sum( plane( myunique(DEG2([11 16 21 22 23]),0) ) ) * plane( DEG1(4) ) * ~( DEG1(4)==DEG2(21) ); 
	TRACK(4) = sum( plane( myunique(DEG2([22 23 24]),0) ) ) * plane( DEG1(5) ) * ~( DEG1(5)==DEG2(23) );
	TRACK(5) = sum( plane( myunique(DEG2([15 20 23 24 25]),0) ) ) * plane( DEG1(6) ) * ~( DEG1(6)==DEG2(25) );
    %
    %
    %%% END 2nd degree neighborhood %%%
    
    
    
    SCORE = SCORE + sum(TRACK);
    if SCORE >= crit
        ISVALID = true;
        return
    else
        ISVALID = false;
    end
    
    
return



function vec = myunique(vec,direct)
    
    dum = vec;idx=1:numel(vec);
    if direct
        [dum, idx] = sort(vec);
    end
    vec = vec( idx( logical(diff(dum)) ) );

return

