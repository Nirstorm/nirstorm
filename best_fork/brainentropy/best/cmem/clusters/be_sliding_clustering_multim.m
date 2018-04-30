function [OPTIONS, SCR, CLS] = be_sliding_clustering_multim(obj, OPTIONS, wtbOPT)
% BE_CLUSTERIZE_MULTIM clusterizes the sources along the data 
% time window in a block fashion. Blocks do not overlap so we make the 
% hypothesis that data is stationary within blocks.  
%
% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.a
%
% INPUTS:
%     - obj                 : MEM obj structure
%     - OPTIONS    			: Structure of parameters (described in be_cmem_solver.m)
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

% Pads data left and right
for ii  = 1 : numel( OPTIONS.mandatory.DataTypes )
    
    if strcmp(OPTIONS.mandatory.pipeline,'cMEM') 
        OPTIONS.automatic.Modality(ii).mspDATA.F = OPTIONS.automatic.Modality(ii).data;
    end
    
    Pad_l   =   OPTIONS.automatic.Modality(ii).mspDATA.F(:,1) * ones( 1 , fix(tW/2) );
    Pad_r   =   OPTIONS.automatic.Modality(ii).mspDATA.F(:,end) * ones( 1 , fix(tW/2) );
    DT{ii}  =   [Pad_l OPTIONS.automatic.Modality(ii).mspDATA.F Pad_r];
end

% Initialize the MSP+clustering
nbS     = size( OPTIONS.automatic.Modality(1).gain,2 );
SCR     = zeros(nbS, tS);
CLS     = zeros(nbS, tS);
initT   = OPTIONS.clustering.MSP_scores_threshold;

% Null hypothesis (for the threshold for the msp scores)
if strcmpi(OPTIONS.clustering.MSP_scores_threshold,'fdr')
    OPTIONS = be_ridges_nullhypothesis(OPTIONS);
end

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
for ii = 1 : tS
    tmV     =   ( 0 : 2*fix(tW/2) ) + ii;
    
    % compute multitrial scores within window ii
    mtSCR = zeros(nbS, 1);
    for jj = 1 : nT
        
        % compute multimodal scores within trial jj
        mmSCR   = zeros(nbS, 1);
        for kk = 1 : numel( OPTIONS.automatic.Modality )
            
             if any(isnan(DT{kk}(:,ii+fix(tW/2))))
                temp = mmSCR;
             else
                [temp, OPTIONS] = be_msp( DT{kk}(:,tmV), OPTIONS.automatic.Modality(kk).gain_struct, OPTIONS);
             end
             
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
    SCR(:,ii)   = mtSCR;
    
    if nargout == 3
        
        if isempty(OPTIONS.optional.clustering.clusters)
            % clusterize dipoles
            OPTIONS.clustering.MSP_scores_threshold     = min(  mean(th) );
            [OPTIONS, temp]                             = be_create_clusters( obj.VertConn, mtSCR, OPTIONS );
            CLS(:,ii)                                   =  temp;
            OPTIONS.clustering.MSP_scores_threshold     = initT;
            
        else
            % either ONE cluster that is repeated (one column in
            % OPTIONS.optional.clustering.clusters) or we thake the ii th
            % column:
            CLS(:,ii)                                   =   OPTIONS.optional.clustering.clusters(:,min( size(OPTIONS.optional.clustering.clusters,2), ii ));
        end
        
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
