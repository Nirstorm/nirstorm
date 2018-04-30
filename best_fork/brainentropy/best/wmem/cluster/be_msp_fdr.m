function [threshold, OPTIONS] = be_msp_fdr(scores, param_Beta, OPTIONS)
%   This function returns the threshold for the msp scores used for 
%   clustering. If OPTIONS.clustering is 'full', the threshold is 0
%   (default); if OPTIONS.clustering is 'fdr', the threshold is found
%   using the fdr criteria based on a baseline.
%   Input:
%   scores : msp scores (obtained from a time windows of interest)
%   param_Beta : parameters of the beta fit of the msp distribution (H0)
%   Output:
%   threshold
%   
% 	Reference
%       JM Lina et al., IEEE 2011.
%   Authors: LATIS team, 2011.
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

%%
if ischar(OPTIONS.clustering.MSP_scores_threshold) && strcmp(OPTIONS.clustering.MSP_scores_threshold,'fdr')
    % Parameters:
    P_VALUE = 0.05;      % p_value for the fdr
    NS = length(scores); % should be equal to the number of dipoles

    % FDR
    p_vH0 = p_valuesH0(scores,param_Beta.alpha,param_Beta.beta);
    scores(:,2)    = p_vH0;
    [p_vout,Isout] = sort(p_vH0);
    scores         = scores(Isout,:);
    scores(:,3)    = (1:NS)/NS*P_VALUE;
    scores_select  = scores(scores(:,3)-scores(:,2)>=0,:);
    if ~isempty(scores_select)
        threshold  = scores_select(end,1);
    else
    % FDR does not distinguish baseline from recordings
    threshold  = 1.0;
    end
elseif isnumeric(OPTIONS.clustering.MSP_scores_threshold)
    threshold = OPTIONS.clustering.MSP_scores_threshold;
else
    threshold = 0.0;
end
end


% ------ functions
function p_v = p_valuesH0 ( s, al, be )
p_v    = 1.0-betainc(max(0,min(1,s)),al,be);
end