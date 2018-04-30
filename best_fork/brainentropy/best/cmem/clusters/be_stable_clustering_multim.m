function [SCR, CLS, OPTIONS] = be_stable_clustering_multim(obj, OPTIONS)
% BE_STABLE_CLUSTERING_MULTIM clusterizes the sources producing
% a spatial clustering stationary in time, while providing
% spatio-temporal evolution of MSP scores.
%
% This method is for data-driven spatio-temporal clustering
% has been described and validated in
% Chowdhury R.A., Lina J.M., Kobayashi E. and Grova C.
% MEG Source Localization of spatially extended generators of epileptic activity:
% Comparing Entropic and Hierarchical Bayesian Approaches.
% PLoS ONE 2013: 8(2):e55969
%
% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.a
%
% INPUTS:
%     - obj        : MEM obj structure
%     - OPTIONS    : Structure of parameters (described in be_main.m)
%        -automatic.Modality    : Structure containing the data and gain matrix for each modality
%          |- data              : Data matrix to process (channels x time samples)
%          |- gain_struct       : Structure returned by be_main_leadfields to give the normalized lead field.
%        -clustering            : Structure containing parcellization parameters
%          |- neighborhood_order: (Set to default 4) Nb. of neighbors used to clusterize
%                                 cortical surface.
%        -optional              : Structure for setting optional parameters
%          |- cortex_vertices   : List of vertices in the cortical mesh.
%        -VertConn              : Neighborhood matrix
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
% Copyright (C) 2011 - MultiFunkIm &  LATIS Team
%
%  Authors: MultiFunkIm team, LATIS team, 2011
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

% Specific options
SVD_threshold = 0.95;

% Needed Data
OPTIONS.clustering.MSP_window = 1;

% test the field ModalityIndex in obj
if ~isfield(obj,'ModalityIndex')
    obj.ModalityIndex = 1;
end

VertConn = obj.VertConn;

% ==== clustering technique
nbS     = size(OPTIONS.automatic.Modality(1, 1).gain,2);        % Nb of sources                                         
tS      = size(OPTIONS.automatic.Modality(1, 1).data,2);        % échantillon temp.                                     
SCR     = zeros(nbS,tS);                                        % nb sources fct of time                               

%% MSP Calculation for each modalities
for ii  = 1 : numel( OPTIONS.mandatory.DataTypes )              %For each modalities
    
    M = OPTIONS.automatic.Modality(ii).data;                    
    Gstruct = OPTIONS.automatic.Modality(ii).gain_struct;       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1- SVD of the DATA MATRIX
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [U,S,V] = svd(M,0);
    s = diag(S);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2- Threshold to identify the signal subspace (min 3 or 95% inertia)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inertia = zeros(1, length(s));                                          % Prealocating inertia to increase the speed
    for i = 1:length(s)
        inertia(i) = sum(s(1:i).^2)./(sum(s.^2));
    end
    
    % cumsum
    q = find(inertia>=SVD_threshold,1);
    
    % Ask for standard display with verbose
    if OPTIONS.optional.verbose
        fprintf('%s, stable clustering: dimension of the signal subspace %d, for inertia > %3.2f\n', OPTIONS.mandatory.pipeline,q,SVD_threshold);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %3- MSP applied on each principal component of the signal subspace
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SPATIO_TEMPORAL EXTENSION of MSP scores in Signal Subspace (Refer to APPENDIX S1 in Chowdhury et al. 2013 Plos One
    
    APM_tot = zeros(nbS,tS);
    
    for i = 1:q
        
        M_data = U(:,i)*S(i,i)*V(:,i)';
        [APM, OPTIONS] = be_msp( M_data, Gstruct, OPTIONS);                 % MSP is calculated here
        
        scale_prob = max(max(abs(V(:,1:q))));
        
        % SPATIO_TEMPORAL VERSION of MSP scores in Signal Subspace
        APM_tot = APM_tot + APM * abs(V(:,i)/scale_prob)';
    end
    
    APM_tot = APM_tot/(max(max(APM_tot)));
    
    SCR = SCR + APM_tot - SCR.*APM_tot;             %To consider every modalities, we should sum these scores (inspired by be_sliding_clustering.)
          
end

% Classification matrix:
% Contains labels ranging from 0 to number of parcels (1 column / time sample) for each sources.

nb_time = size(SCR,2);
CLS     = zeros(nbS, nb_time);

for i = 1:nb_time
    [OPTIONS, CLS(:,i),cellstruct_cls] = be_create_clusters(VertConn, mean(SCR,2), OPTIONS );
end

end