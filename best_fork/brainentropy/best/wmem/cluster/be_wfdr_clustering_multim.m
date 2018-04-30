function [CLS, SCR, OPTIONS] = be_wfdr_clustering_multim(obj, OPTIONS)
% BE_FDR_CLUSTERIZE_MULTIM clusterizes the sources along the data 
% time window in a wavelet representation. Only tf-boxes with global energy with
% more than 0.1% of the variance (at each scale) is kept in the model.  
%
% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.a
%
% INPUTS:
%     - WM      : Data matrix (wavelet rep.) to be reconstructed using MEM
%     - Gstruct : Structure returned by normalize_gain. See normalize_gain
%                 for more info.
%     - Nm      : Neighborhood matrix
%     - OPTIONS : Structure of parameters (described in be_memsolver_multiM.m)
%
% OUTPUTS:
%     - Results : Structure
%          |- OPTIONS       : Keep track of parameters
%          |- CLS           : Classification matrix. Contains labels
%          |                  ranging from 0 to number of parcels (1 column
%          |                  by time sample) for each sources.
%          |- SCR           : MSP scores matrix (same dimensions a CLS).
%
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

    % H. Mahkoum, JM Lina's appraoch:
    if OPTIONS.optional.verbose
        fprintf('%s, fdr clustering ...\n', OPTIONS.mandatory.pipeline);
    end

% Initialize the MSP+clustering
nbb = size(OPTIONS.automatic.selected_samples,2);
nbS = size(OPTIONS.automatic.Modality(1).gain,2);
% SCR contains the msp scores for all sources (lines), in each t-f box
% (columns):
SCR = zeros(nbS, nbb);
CLS = zeros(nbS, nbb);
% there is one MSP threshold per t-f box: 
THR = ones(1,nbb);
Mod_in_box  = zeros(1,nbb);
initTH      = OPTIONS.clustering.MSP_scores_threshold;

if OPTIONS.optional.verbose
    fprintf('%s, MSP-based clustering -> boxes(threshold,nb. clusters) follow:\n\t', OPTIONS.mandatory.pipeline); 
end 
% Loops on all data sub-windows and modalities 
for ii = 1 : nbb
    % the box may belong to one of the modalities or both.
    % if 0 the box_in_Mod will not consider the modality 
    % if 1 compute multimodal scores within the box ii
    combSCR = zeros(nbS, 1);
    th      = [1 1];
    for jj = 1 : numel(OPTIONS.mandatory.DataTypes)
        box_in_Mod = ii;% OPTIONS.automatic.selected_samples(1,ii);
        if box_in_Mod
%           ST = OPTIONS.automatic.Modality(jj).selected_jk(4,box_in_Mod);  % previously 5 instead of 4
%           ND = OPTIONS.automatic.Modality(jj).selected_jk(5,box_in_Mod);  % previously 6 instead of 5         
            ST = OPTIONS.automatic.selected_samples(4,box_in_Mod);  % previously 5 instead of 4
            ND = OPTIONS.automatic.selected_samples(5,box_in_Mod);  % previously 6 instead of 5  
            [scores, OPTIONS] = be_msp(OPTIONS.automatic.Modality(jj).data(:,ST:ND), OPTIONS.automatic.Modality(jj).gain_struct, OPTIONS);
            
            if isempty(OPTIONS.optional.clustering.clusters)
                % MSP scores threshold
                if strcmp(OPTIONS.clustering.MSP_scores_threshold, 'fdr')
                    [th(jj), OPTIONS] = be_msp_fdr(scores, OPTIONS.automatic.Modality(jj).paramH0, OPTIONS);
                elseif isnumeric(OPTIONS.clustering.MSP_scores_threshold)
                    th(jj) = OPTIONS.clustering.MSP_scores_threshold;
                end
            end
            
            combSCR = combSCR + scores - combSCR.*scores;
            Mod_in_box(1,ii) = Mod_in_box(1,ii)+jj;
        end
    end
    SCR(:,ii) = combSCR;
    THR(1,ii) = min(th);    
    
    if isempty(OPTIONS.optional.clustering.clusters)
        % clusterize dipoles
        OPTIONS.automatic.Mod_in_boxes          = [Mod_in_box ; THR];
        OPTIONS.clustering.MSP_scores_threshold = THR(1,ii); % to force the threshold
        [OPTIONS, temp]                         = be_create_clusters( obj.VertConn, SCR(:,ii), OPTIONS );
        OPTIONS.clustering.MSP_scores_threshold = initTH; % to restore the fdr option
        % the cluster config of the sources
        CLS(:,ii) =  temp;
        if OPTIONS.optional.verbose
            fprintf('%3d(%.2f,%3d),', OPTIONS.automatic.selected_samples(1,ii), THR(1,ii), max(temp)); 
            if mod(ii,5)==0, fprintf('\n\t'), end
        end
    else
        %CLS(:,ii) = OPTIONS.optional.clustering(:, min(size(OPTIONS.optional.clustering,2),ii) );
        CLS(:,ii) = OPTIONS.optional.clustering.clusters(:, min(size(OPTIONS.optional.clustering.clusters,2),ii) );
    end
        
end
if OPTIONS.optional.verbose
    fprintf('\n%s, %d boxes have been selected.\n',OPTIONS.mandatory.pipeline,nbb); 
end

return