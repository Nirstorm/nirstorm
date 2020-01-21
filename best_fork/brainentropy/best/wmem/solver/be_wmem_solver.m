function [Results, OPTIONS] = be_wmem_solver(HeadModel, OPTIONS, Results)
% MEMSOLVER: Maximum Entropy on the Mean solution.
%
% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.a
%
% INPUTS:
%     - HeadModel  : Brainstorm head model structure
%     - OPTIONS    : Structure of parameters (described in be_sourceimaging.m)
%          |- spatial_smoothing : optional spatial constraint on the
%                                 variance of the active state function, 
%                                 for each cluster. Set to 1 to activate.
%          |- MSP_R2_threshold  : threshold (between 0 and 1) to filter 
%                                 data and forward operator in the msp 
%                                 calculation
%          |- Modality          : Cell array of strings that identifies which
%                                 modalities are selected
%          |- MSP_window        : number of samples to use for calculating
%                                 the MSP. The window should be small 
%                                 enough to satisfy assumption of 
%                                 stationarity within that window.
%          |- stable_clusters   : Set to 1 to used C. Grova's clustering
%                                 method
%          |- NoiseCov          : Noise variance-covariance matrix. Should
%                                 be available from Brainstorm or can be 
%                                 estimated and provided separately
%          |- MSP_scores_threshold :Threshold of the MSP scores (between 0 
%                                   and 1). Set to 0 to keep all the scores.
%          |- active_var_mult   : multuplier for the variance-covariance of
%                                 the active state function
%          |- inactive_var_mult : multuplier for the variance-covariance of
%                                 the inactive state function (default: 0
%                                 i.e. a dirac)
%          |-alpha_threshold    : Threshold on the cluster probabilities of
%                                 being active. Clusters whose probas are
%                                 lower than this threshold will be turned
%                                 off
%          |- Optim_method      : Choice of optimization routine. By 
%                                 default uses minFunc (Mark Schmidt's
%                                 implementation). Specify 'fminunc' to use
%                                 Matlab optimization toolbox routine.
%          |- alpha             : Option to specify the alpha probabilies 
%                                 manually
%          |- Data              : Data to process (channels x time samples)
%          |- Units             : Structure for rescaling units
%          |- neighborhood_order: Nb. of neighbors used to cluterize 
%                                 cortical surface
%          |- alpha_method   :  Selection of the method to initialize the
%                               probability (alpha) of activation of the 
%                               clusters  (alpha).
%          |                   1 = average of the scores of sources in a cluster 
%          |                   2 = maximum of the scores of sources in a cluster 
%          |                   3 = median of the scores of sources in a cluster 
%          |                   4 = equiprobable activity (0.5)
%          |- active_mean_method:  Selection of the method to initialize 
%                                  the mean of the "active state" law (mu).
%          |                   1 = null hypothesis
%          |                   2 = Null activity (0)
%
% OUTPUTS:
%     - Results : Structure
%          |- ImageGridAmp  : Source activity (sources x times)
%          |- ImagingKernel : Not computed ([])
%
% ==============================================
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2011
%  ref [1]. Lina and al., wavelet-based localization of oscillatory sources
%  from MEG data, IEEE TBME (2012)
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
%%
global MEMglobal

% pipeline starts here:
if OPTIONS.optional.verbose
    fprintf('\n\n===== pipeline wMEM\n');
end        
        
%% Useful variables
obj.ImageGridAmp = []; 
if ~isfield(HeadModel, 'vertex_connectivity')
    [OPTIONS, obj.VertConn] = be_vertex_connectivity(HeadModel, OPTIONS);
else
    obj.VertConn = HeadModel.vertex_connectivity;
end

if isempty(OPTIONS.optional.clustering) && ( isempty(obj.VertConn) || diff(size(obj.VertConn)) ) % CETTE PARTIE N'EST PAS CORRECTE: || size(obj.VertConn,1)~=size(HeadModel.Gain,2)
    fprintf('wMEM solver error : no vertex con.nectivity matrix available.\n');
    return
end


%% ===== Comment ===== %%
OPTIONS.automatic.Comment       =   OPTIONS.optional.Comment;
if strcmp( OPTIONS.automatic.Comment(1:3), 'MEM' )
    OPTIONS.automatic.Comment   =   ['w' OPTIONS.optional.Comment];
end

%% ===== Channels ===== %% 
% we retrieve the channels name and the data
[OPTIONS, obj]  = be_main_channel(HeadModel, obj, OPTIONS);

%% ===== Sources ===== %% 
% we verify that all sources in the model have good leadfields
[OPTIONS, obj]  = be_main_sources(obj, OPTIONS);

%% ===== Pre-whitening of the data ==== %%
% it uses empty-room data if available
% if PlOS one : nothing is done here
% [OPTIONS] = be_prewhite(OPTIONS);

%% ===== Pre-process the leadfield(s) ==== %% 
% we keep leadfields of interest; we compute svd of normalized leadfields
[OPTIONS, obj] = be_main_leadfields(obj, OPTIONS);

%% ===== Normalization ==== %% 
% we absorb units (pT, nA) in the data, leadfields; we normalize the data
% and the leadfields
[OPTIONS] = be_normalize_and_units(OPTIONS);

%% ===== Null hypothesis (for the threshold for the msp scores)
% from the baseline, compute the distribution of the msp scores. 
% More details in Appendix B in ref[1]
OPTIONS = be_model_of_null_hypothesis(OPTIONS);

%% ===== Data Processing (and noise) ===== %%
% for the data: normalization/wavelet/denoise
[OPTIONS, obj] = be_wdata_preprocessing(obj, OPTIONS);
if OPTIONS.optional.display
    [obj.hfig obj.hfigtab] = be_display_time_scale_boxes(obj,OPTIONS);
end

%% ===== Double to single precision  ===== %%
[OPTIONS] = be_switch_precision(OPTIONS, 'single');

%% ===== Clusterize cortical surface ===== %%
% from the msp scores, clustering of the cortical mesh:
[OPTIONS, obj] = be_main_clustering(obj, OPTIONS);

%% ===== Single to double precision  ===== %%
[OPTIONS] = be_switch_precision(OPTIONS, 'double');

%% ===== pre-processing for spatial smoothing (Green mat. square) ===== %%
% matrix W'W from the Henson paper
[OPTIONS, obj.GreenM2] = be_spatial_priorw( OPTIONS, obj.VertConn);

%% ===== Solve the MEM ===== %%
% the main inverse process
[obj, OPTIONS] = be_main_wmem(obj, OPTIONS);
if OPTIONS.optional.display
    be_display_entropy_drops(obj,OPTIONS);
end

% Results (full temporal sequence)
Results = struct(...
    'ImageGridAmp',     obj.ImageGridAmp, ...
    'ImagingKernel',    [], ...
    'MEMoptions',       OPTIONS, ...
    'MEMdata',          obj);

return


