function [groups, ctf] = be_group_ridges(ridges, duration, energy)
%BE_GROUP_RIDGES identifies connected ridge lines from a binary mask
%
%   INPUT:
%       -   ridges  :   local maxima mask
%       -   duration:   duration threshold on ridge lines
%       -   energy  :   energy threshold on ridge lines
%
%   OUTPUTS:
%       -   groups  :   supra-thresholds ridge lines
%       -   ctf     :   average frequencies of ridges lines
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

%
%% Algorithme glouton
candidats = find(ridges > 0);

if ~numel(candidats)
    groups  =   {};
    ctf     =   {};
else
    NRJ     =   sum( ridges(:).^2 );
    aslice  =   zeros(size(ridges));
    aslice(candidats)   =   1;
    if duration>0
        BW  =   bwareaopen(aslice,duration);
        BW  =   bwlabel(BW);
    else
        BW  =   bwlabel(aslice);
    end
    S       =   regionprops(BW,ridges,'Area','PixelValues','PixelIdxList');
    ctf     =   cellfun(@(x) mean(mod(x,size(ridges,1))),{S.PixelIdxList},'uni',false);
    lgts    =   [S.Area];
    nrjs    =   cell2mat(cellfun(@(a) sum(a.^2),{S.PixelValues},'uni',false));
    strs    =   cell2mat( cellfun(@(a) sum(a),{S.PixelValues},'uni',false) );
    groups  =   {S.PixelIdxList};
    
    if numel(lgts)
        crit        =   lgts .* strs;
        [dum, indx] =   sort(crit);
        groups(:)   =   groups(fliplr(indx));
        ctf(:)      =   ctf(fliplr(indx));
        
        if numel(nrjs)>1&&sum(nrjs> energy*NRJ);
            [dum, II] = find(cumsum(nrjs)>energy*NRJ);
            if numel(II)>1
                groups( II(2) : end )   = [];
                ctf( II(2) : end )      = [];
            end
        end
    end
    
    groups  =   cellfun(@(a) a', groups, 'uni', 0);
    
end

return