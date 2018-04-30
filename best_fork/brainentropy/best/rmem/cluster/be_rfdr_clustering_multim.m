function [CLS, SCR, OPTIONS] = be_rfdr_clustering_multim(obj, OPTIONS)
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
%     - OPTIONS    : Structure of parameters (described in be_memsolver_multiM.m)
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

    % Zerouali's appraoch:
    if OPTIONS.optional.verbose
        fprintf('%s, fdr clustering ...\n', OPTIONS.mandatory.pipeline);
    end
% Load raw data
T	=   OPTIONS.automatic.Modality(1).mspDATA.Time;
Frq =   OPTIONS.automatic.Modality(1).mspDATA.FRQs;

% Initialize the MSP+clustering
nbb     = size(OPTIONS.automatic.Modality(1).data, 2);
nbS     = size(OPTIONS.automatic.Modality(1).gain, 2);
nMod    = numel( OPTIONS.mandatory.DataTypes );
nbL     = numel(Frq);
SCR     = zeros(nbS, nbb);
CLS     = zeros(nbS, nbb);
THR     = ones(1, nMod);
%rdg_in_Mod = [1 1];
%Mod_in_rdg = zeros(1,nbb);
init       = OPTIONS.clustering.MSP_scores_threshold;

% Get ridges central frequency and adapt fdr window
%% ===== Null hypothesis (for the threshold for the msp scores)
OPTIONS     =   be_ridges_nullhypothesis(OPTIONS);

fprintf('%s-MSP-based clustering -> ridges(threshold,nb. clusters) follow:\n\t', OPTIONS.mandatory.pipeline); 
 

                
% Loops on all data sub-windows and modalities 
for ii = 1 : nbb
    % the box may belong to one of the modalities or both
    % the box_in_Mod will follow the progression
    % compute multimodal scores within the box ii
    rdgID   = OPTIONS.automatic.selected_samples(ii);
    rdgFR   = Frq( mod( rdgID, nbL ) );
    if strcmp(OPTIONS.clustering.MSP_scores_threshold, 'fdr') 
        winW    = fix( OPTIONS.ridges.cycles_in_window / rdgFR * OPTIONS.automatic.sampling_rate ) + 1;    
    elseif isnumeric(OPTIONS.clustering.MSP_scores_threshold) 
        winW    = OPTIONS.clustering.MSP_window;  
    end
    winW    = winW - mod(winW,2);
    combSCR = zeros(nbS, 1);
    th      = ones( 1, nMod );
    for jj = 1 : nMod
        if OPTIONS.ridges.fdrRdgMod{jj}(rdgID)
            ST = max( ii-winW/2, 1 );
            ND = min( ii+winW/2, size(OPTIONS.automatic.Modality(jj).mspDATA.F,2) );       
            [scores, OPTIONS] = be_msp( OPTIONS.automatic.Modality(jj).mspDATA.F(:,ST:ND), ...
                                OPTIONS.automatic.Modality(jj).gain_struct, OPTIONS);
            
            if strcmp(OPTIONS.clustering.MSP_scores_threshold, 'fdr') 
                [th(jj), OPTIONS] = be_msp_fdr(scores, OPTIONS.automatic.Modality(jj).paramH0, OPTIONS);
            elseif isnumeric(OPTIONS.clustering.MSP_scores_threshold) 
                th(jj) = OPTIONS.clustering.MSP_scores_threshold;
            end
            combSCR = combSCR + scores - combSCR.*scores;
%             rdg_in_Mod(1,jj) = rdg_in_Mod(1,jj)+1;
%             Mod_in_rdg(1,ii) = Mod_in_rdg(1,ii)+jj;
        end
    end
    SCR(:,ii) = combSCR;
    THR(1,ii) = min(th);    

    if isempty(OPTIONS.optional.clustering.clusters)
        % clusterize dipoles
        % OPTIONS.automatic.Mod_in_ridges = [Mod_in_rdg ; THR];
        OPTIONS.clustering.MSP_scores_threshold = THR(1,ii); % to force the threshold
        [OPTIONS, temp] = be_create_clusters( obj.VertConn, SCR(:,ii), OPTIONS );
        OPTIONS.clustering.MSP_scores_threshold = init; % to restore the fdr option
        % the cluster config of the sources
        CLS(:,ii) =  temp;
        fprintf('%3d(%.2f,%3d),', OPTIONS.automatic.selected_samples(1,ii), THR(1,ii), max(temp)); 
        if mod(ii,7)==0, fprintf('\n\t'), end
    else
        CLS(:,ii) = OPTIONS.optional.clustering.clusters(:, min( size(OPTIONS.optional.clustering.clusters,2), ii) );
    end
    
end
    fprintf('\n%d boxes have been selected.\n',nbb); 
    OPTIONS.automatic.selected_samples = [];
return