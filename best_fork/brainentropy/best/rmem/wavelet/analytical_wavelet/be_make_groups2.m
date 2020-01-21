function [GROUPS, CTF, OPTIONS] = be_make_groups2( RIDGES, OPTIONS, varargin )
% BE_MAKE_GROUPS2 groups ridge points into ridge lines 
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


ID = 1;
if numel(varargin) == 1 && isempty(varargin{1});
    ID = varargin{1};
end

[RIDGES, params, OPTIONS] = be_focus_on_window( RIDGES, OPTIONS, 1, ID );

% Cleans the ridges plan - local maxima of the ridge map at each sample
if OPTIONS.ridges.skim_map
    mask    = be_skim_ridges(RIDGES);
    RIDGES  = RIDGES .* mask;
end

% Make the groups
[GROUPS, CTF]   = be_group_ridges(RIDGES, OPTIONS.ridges.min_duration, ...
                  OPTIONS.ridges.energy_threshold);
GROUPS          = cellfun(@(a) ceil(a/params.nbL2+params.smp1-2)* ...
                  params.nbL+params.initF2+ mod(a-1, params.nbL2)+1, ...
                  GROUPS, 'uni', false);
GROUPS          = cellfun(@(a) [a(:)'], GROUPS, 'uni', false);
CTF             = cellfun(@(a) ceil(a/params.nbL2+params.smp1-2)* ...
                  params.nbL+params.initF2+ mod(a-1, params.nbL2)+1, ...
                  CTF, 'uni', false);
CTF             = cellfun(@(a) mod(a-1, params.nbL)+1, CTF, 'uni', false);

return







