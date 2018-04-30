function [OPTIONS] = be_avg_reference(OPTIONS)
% This function average reference the data of every channels,
% along every time samples
%
%% ==============================================
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2015
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

for ii = 1:numel(OPTIONS.mandatory.DataTypes)
         
    % Compute mean along channels for every time sample
    mu = mean(OPTIONS.automatic.Modality(ii).data);    
    muM = ones(size(OPTIONS.automatic.Modality(ii).data,1),1)*mean(OPTIONS.automatic.Modality(ii).data);

    % subtract from data
    OPTIONS.automatic.Modality(ii).data = OPTIONS.automatic.Modality(ii).data - muM;
    
end
return