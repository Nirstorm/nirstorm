function  display_entropy_drops(obj,OPTIONS)
%DISPLAY_ENTROPY_DROPS displays a graphical window with MEM entropy drops at 
% each  iteration
%
%   INPUTS:
%       -   obj
%       -   OPTIONS
%
%   OUTPUTS:
%      
%% ==============================================   
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2011
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

warning('off','all');
    onglet = uitab(obj.hfigtab,'title','ENTROPY');
    
    hpc = uipanel('Parent', onglet, ...
                  'Units', 'Normalized', ...
                  'Position', [0.01 0.01 0.98 0.98], ...
                  'FontWeight','demi');

    set(hpc,'Title',[' Entropy drops'],'FontSize',8);

    ax = axes('parent',hpc, ...
    'outerPosition',[0.01 0.01 0.98 0.98]);
   
    plot(ax,log(abs(OPTIONS.automatic.entropy_drops)),'-k','linewidth',2);
    xlabel(ax,'wavelet coefficient (descreasing energy)'); 
    ylabel(ax,'log(|entropy|)'); 
return