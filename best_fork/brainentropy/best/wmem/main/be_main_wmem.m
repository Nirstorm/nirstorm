function [obj, OPTIONS] = be_main_wmem(obj, OPTIONS)
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

	if strcmp(OPTIONS.mandatory.pipeline,'wMEM')
        	obj = be_fusion_of_modalities(obj.data,obj,OPTIONS);
        	[obj.ImageGridAmp, OPTIONS] = be_launch_mem(obj, OPTIONS);
            
            % si j=0, on remplace obj.data = obj.scaling_data, on corrige ? la matrice de variance covariance
            % et relancer le be_launch_mem
            % OPTIONS_scaling = OPTIONS. 
            % OPTIONS_scaling.automatic.selected_samples = selection qui correspond au coeff d'echelle;
            % [obj.scaling, OPTIONS_scaling] = be_launch_mem(obj, OPTIONS_scaling);
            % obj.ImageGridAmp = obj.ImageGridAmp + obj.scaling;
            
        	if ~isempty(obj.ImageGridAmp)
            	obj.ImageGridAmp = obj.ImageGridAmp(:,obj.info_extension.start:obj.info_extension.end);
            end
        	OPTIONS.automatic.Comment = [OPTIONS.automatic.Comment ' DWT(j' num2str(OPTIONS.wavelet.selected_scales) ')'];
	else
		disp('ERROR in calling wMEM'); return
	end

end
