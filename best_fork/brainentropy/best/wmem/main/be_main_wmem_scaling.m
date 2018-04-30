function [obj_scaling] = be_main_wmem_scaling(obj, OPTIONS)
% BE_MAIN_MEM sets the appropriate options for the MEM 
% accroding to the chosen MEM pipeline
%
%   INPUTS:
%       -   obj
%       -   OPTIONS
%
%   OUTPUTS:
%       -   OPTIONS
%       - obj
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


for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
    data = zeros(size(obj.data{ii}));
    data(:,1:size(obj.scaling_data{ii},2)) = obj.scaling_data{ii};
    obj.data{ii} = data;
end

Nj = fix(log2(size(data,2)));
max_scale = OPTIONS.automatic.scales(1,end);
iBoxes_max_scale = OPTIONS.automatic.selected_samples(2,:)==max_scale;
selected_k = OPTIONS.automatic.selected_samples(3,iBoxes_max_scale);

selected_samples = zeros(6,sum(iBoxes_max_scale)); % set to 6 because of line 49 in be_launch_mem; 3 is enough?
selected_samples(1,:) = selected_k+1;
selected_samples(2,:) = Nj;
selected_samples(3,:) = selected_k;

OPTIONS.automatic.selected_samples = selected_samples;
[obj_scaling, OPTIONS_scaling] = be_launch_mem(obj, OPTIONS);
% we do not keep this OPTIONS in the output
% WE HAVE TO THINK ABOUT THE CLUSTERING ASSOCIATED TO THIS PART OF THE
% ESTIMATOR


end
