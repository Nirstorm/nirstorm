function [OPTIONS, SCR, CLS] = be_stepwise_clustering_multim(obj, OPTIONS, wtbOPT)
% BE_STEPWISE_CLUSTERING_MULTIM clusterizes the sources along the data 
% time window in a block fashion. Blocks do not overlap so we make the 
% hypothesis that data is stationary within blocks.  
%
% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.
%
% INPUTS:
%     - obj                 : MEM obj structure
%     - OPTIONS             : Structure of parameters (described in be_cmem_solver.m)
%     - wtbOPT              : obsolete...
%
% OUTPUTS:
%     - Results : Structure
%          |- OPTIONS       : Keep track of parameters
%          |- SCR           : MSP scores matrix (same dimensions as CLS).
%          |- CLS           : Classification matrix. Contains labels
%          |                  ranging from 0 to number of parcels (1 column
%          |                  by time sample) for each sources.
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


tW = OPTIONS.clustering.MSP_window; % window size
[dum, tS, nT] = size( OPTIONS.automatic.Modality(1).data ); % data length

% Creates pointers for data segmentation into sub-windows of length tW
if tS <= tW
    vec = [0 tS];
else
    divider = tS/tW;
    if divider==fix(divider)
        vec = [0 [1:divider] * tW];
    else
        vec = [0 (1:fix(divider))*tW tS];
    end
end

% Initialize the MSP+clustering
nbS     = size( OPTIONS.automatic.Modality(1).gain,2 );
SCR     = zeros(nbS, tS);
CLS     = zeros(nbS, tS);
initT   = OPTIONS.clustering.MSP_scores_threshold;

% Null hypothesis (for the threshold for the msp scores)
OPTIONS = be_ridges_nullhypothesis(OPTIONS);

% Progress bar
hmem = [];
if ~isempty(wtbOPT)
    hmem    = wtbOPT(1);
    st      = wtbOPT(2);
    dr      = wtbOPT(3);
    prg1    = 0;
    prg2    = 0;
end

% Loops on all data sub-windows and modalities 
for ii = 2 : numel(vec)
    ST = vec(ii-1)+1;
    ND = vec(ii);
    
    % compute multitrial scores within window ii
    mtSCR = zeros(nbS, 1);
    for jj = 1 : nT
        
        % compute multimodal scores within trial jj
        mmSCR   = zeros(nbS, 1);
        for kk = 1 : numel( OPTIONS.automatic.Modality )
            [temp, OPTIONS] = be_msp( OPTIONS.automatic.Modality(kk).data(:,ST:ND), OPTIONS.automatic.Modality(kk).gain_struct, OPTIONS);
            
            if ischar(OPTIONS.clustering.MSP_scores_threshold) && strcmpi(OPTIONS.clustering.MSP_scores_threshold, 'fdr')
                if isempty(OPTIONS.automatic.Modality(kk).paramH0)
                    %% ===== Null hypothesis (for the threshold for the msp scores)
                    windS   = OPTIONS.clustering.MSP_window;
                    OPTIONS.temporary.FDR_max_window_size = fix( 1.2 * windS ) + 1;
                    OPTIONS.temporary.FDR_min_window_size = fix( 0.8 * windS ) + 1;
                    OPTIONS = be_model_of_null_hypothesis(OPTIONS);
                end
                [th(jj,kk), OPTIONS] = be_msp_fdr(temp, OPTIONS.automatic.Modality(kk).paramH0(jj), OPTIONS);
                
            else
                th = initT;
            end
            mmSCR = mmSCR + temp - mmSCR.*temp;
            
            % update progress bar
            if hmem
                prg  = ( (jj - 1  + kk/numel( OPTIONS.automatic.Modality )) / nT ) / (numel(vec)-1);
                prg1 = prg2 + round( (st + dr * prg/2) * 100 );
                waitbar(prg1/100, hmem, ['Step 1/2 : Running cortex parcellization ... ' num2str(prg1) ' % done']);
            end
        end
    
        mtSCR = mtSCR + mmSCR / nT;
         
    end    
    SCR(:,ST:ND)                    = mtSCR * ones(1, ND-ST+1);
    
    if nargout == 3
        % clusterize dipoles
        OPTIONS.clustering.MSP_scores_threshold     = min(  mean(th) );
        [OPTIONS, temp]                             = be_create_clusters( obj.VertConn, SCR(:,ST), OPTIONS );
        CLS(:,ST:ND)                                =  temp * ones(1, ND-ST+1);
        OPTIONS.clustering.MSP_scores_threshold     = initT;
         % update progress bar
         if hmem
             prg2 = prg1 + round( (st + dr * prg/2) * 100 );
             waitbar(prg2/100, hmem, ['Step 1/2 : Running cortex parcellization ... ' num2str(prg2) ' % done']);
         end
    end
end


if isfield(OPTIONS, 'temporary')
    OPTIONS = rmfield(OPTIONS, 'temporary');
end

return