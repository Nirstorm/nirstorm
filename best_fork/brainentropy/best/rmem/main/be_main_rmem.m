function [obj, OPTIONS] = be_main_MEM(obj, OPTIONS)
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


	if strcmp(OPTIONS.mandatory.pipeline,'rMEM')
        	new_obj                     =   be_fusion_of_modalities( [], obj, OPTIONS);
            new_obj.noise_var           =   real(new_obj.noise_var);
        	[obj.ImageGridAmp, OPTIONS] =   be_launch_mem(new_obj, OPTIONS); 
        	% - imaginary part 
        	if ~isempty(new_obj.idata)
            	new_obj                 =   be_fusion_of_modalities( [], obj, OPTIONS);
                if ~isreal( new_obj.noise_var )
                    new_obj.noise_var  	=   imag(new_obj.noise_var);
                end
                new_obj.data            =   new_obj.idata;
            	obj.ImageGridAmp        =   obj.ImageGridAmp + 1i * be_launch_mem(new_obj, OPTIONS);
        	end
	else
		disp('ERROR in calling rMEM'); return
	end

end
