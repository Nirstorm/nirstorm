function [groups, CtF] = be_make_groups( ridges, time, OPTIONS )
% MAKE_GROUPS groups ridge points into ridge lines 
%
% Auteur: LATIS
% Date: Feb 2, 2011 
%
% DESCRIPTION
%   Cette fonction est un algorithme glouton qui explore un plan de 
%   ridge superpos�s et suit les lignes de ridges pour former des crete 
%   qui sont group�es dans la structure obj.groups. Les groupes sont 
%   class�s par taille
% 
%   INPUTS:
%       -   ridges  :   les cartes de ridges de chaque capteur
%       -   time    :   time vector
%       -   OPTIONS :   options structure
%
%   OUTPUTS:
%       -   groups  :   structure qui stocke les indices des ridges qui
%                       forment chacune des lignes
%       -   CtF     :   average frequency of each ridge line
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

% Get appropriate time segment
[nbL, nbS]  = size( ridges );
frqA        = OPTIONS.wavelet_freqs_analyzed;
[dum, smp1] = min( abs(OPTIONS.time_analyzed(1) - time) );
[dum, smp2] = min( abs(OPTIONS.time_analyzed(end) - time) );
% Get appropriate frequency range
[dum, frq1] = min( abs(OPTIONS.ridges_frequency_range(1) - frqA) );
[dum, frq2] = min( abs(OPTIONS.ridges_frequency_range(2) - frqA) );

ridges( :, [1:max(1,smp1-1) min(smp2+1,nbS):end] ) = 0; 
ridges( [1:max(1,frq2-1) min(frq1+1,nbL):end], : ) = 0;

ridges( ridges < OPTIONS.ridges_strength_threshold ) = 0;

if OPTIONS.ridges_skim_map
    % Cleans the ridges plan - local maxima of the ridge map at each sample
    mask = be_skim_ridges(ridges);
    ridges = ridges .* mask;
end

NRJ = sum( ridges(:).^2 );


%% Algorithme glouton
candidats = find(ridges > 0);

if ~numel(candidats)
    groups = {};
    CtF = {};
else
    aslice=zeros(size(ridges));
    aslice(candidats)=1;
    if OPTIONS.ridges_duration_threshold>0
        BW=bwareaopen(aslice,OPTIONS.ridges_duration_threshold);
        BW=bwlabel(BW);
    else
        BW=bwlabel(aslice);
    end
    S=regionprops(BW,ridges,'Area','PixelValues','PixelIdxList');
    CtF=cellfun(@(x) mean(mod(x,nbL)),{S.PixelIdxList},'uni',false);
    lgts=[S.Area];
    nrjs=cell2mat(cellfun(@(a) sum(a.^2),{S.PixelValues},'uni',false));
    strs = cell2mat( cellfun(@(a) sum(a),{S.PixelValues},'uni',false) );
    
    groups={S.PixelIdxList};
    
    if numel(lgts)
        crit = lgts .* strs;
        [dum, indx] = sort(crit);
        groups(:) = groups(fliplr(indx));
        CtF(:) = CtF(fliplr(indx));
        
        if numel(nrjs)>1&&sum(nrjs> OPTIONS.ridges_energy_threshold*NRJ)
            [dum, II] = find(cumsum(nrjs)>OPTIONS.ridges_energy_threshold*NRJ);
            if numel(II)>1
                groups( II(2) : end ) = [];
                CtF( II(2) : end ) = [];
            end
        end
    end
end



