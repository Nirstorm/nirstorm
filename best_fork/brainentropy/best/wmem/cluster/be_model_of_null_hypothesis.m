function OPTIONS = be_model_of_null_hypothesis(OPTIONS)
% BE_MODEL_OF_NULL_HYPOTHESIS models the distribution of the MSP scores 
% under the null hypothesis. A baseline is mandatory for this function and
% must be inputted approrpriately in the OPTIONS.
% 
%   OPTIONS = be_model_of_null_hypothesis(OPTIONS) returns the parameters of the       
%   beta law that best models the distribution of MSP scores of the baseline. 
%   
%   INPUTS:
%       -   OPTIONS :   structure containig all BEst options
%
%   OUTPUTS:
%       -   OPTIONS :   structure with new fields (clustering.MSP_H0)
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

nMod = numel(OPTIONS.mandatory.DataTypes);

is_threshold = ischar(OPTIONS.clustering.MSP_scores_threshold);
is_fdr = strcmp(OPTIONS.clustering.MSP_scores_threshold,'fdr');

if (is_threshold && is_fdr && isempty(OPTIONS.optional.clustering.clusters))
    if OPTIONS.optional.verbose
        fprintf('%s, fdr criteria for clustering\n', OPTIONS.mandatory.pipeline);
    end
    
    for ii = 1 : nMod
        if ~isempty(OPTIONS.automatic.Modality(ii).baseline)
        BSL = OPTIONS.automatic.Modality(ii).baseline*OPTIONS.automatic.Modality(ii).units.Data_units; 
        [param, OPTIONS] = be_msp_h0beta_fit(BSL, OPTIONS.automatic.Modality(ii).gain_struct, OPTIONS);
        OPTIONS.automatic.Modality(ii).paramH0.alpha = param.alpha;
        OPTIONS.automatic.Modality(ii).paramH0.beta  = param.beta;
        if OPTIONS.optional.verbose, fprintf('\t%s : %5.2f (alpha), %5.2f (beta)\n', OPTIONS.automatic.Modality(ii).name,param.alpha,param.beta); end;
        else
        if OPTIONS.optional.verbose, fprintf('\tNo fdr possible'); end;
        end
    end
    
end

end