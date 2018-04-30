function [slice] = be_skim_ridges(slice)
%BE_SKIM_RIDGES filters superimposed ridge planes in order to extract ridge 
% lines
%
%   INPUTS:
%       -   slice   :   superimposed ridge plane
%
%   OUTPUTS:
%       -   slice   :   skimmed plane
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

    cands = find(slice>0);
    while ~isempty(cands)
        CND = cands(1);
        [Dwn] = scan_down(CND, slice);
        vec = [CND Dwn];
        [dum, ind] = max( slice(vec) );
        slice(vec) = 0;
        slice( vec(ind) ) = 1;
        
        [dum, ind] = ismember(vec, cands);
        cands(ind) = [];
    end

return



function [Dwn] = scan_down(CND, slice)

    stop = 0;
    
    [M,N] = size(slice);
     [m,n] = ind2sub([M,N], CND); 
    aft = sort( (n-1)*M+m+1:n*M ); 
    
    % SCAN DOWN
    stop = 0;
    Dwn = [];
    while ~stop && ~isempty(aft)
        if slice( aft(1) )>0
            Dwn = [Dwn aft(1)];
            aft(1) = [];
        else
            stop = 1;
        end
    end    

return