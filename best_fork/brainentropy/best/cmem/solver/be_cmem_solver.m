function [Results, OPTIONS] = be_cmem_solver(HeadModel, OPTIONS, Results)
% cMEMSOLVER: Maximum Entropy on the Mean solution.
%
% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.a
%
% INPUTS:
%     - HeadModel  : Brainstorm head model structure
%     - OPTIONS    : Structure of parameters (described in be_main.m)
%          -solver          : Parameters for MEM solver
%               |- spatial_smoothing : optional spatial constraint on the
%                                 variance of the active state function, 
%                                 for each cluster. Set to 1 to activate.
%               |- NoiseCov          : Noise variance-covariance matrix. 
%                                      Should be available from Brainstorm or can be 
%                                      estimated and provided separately
%               |- active_var_mult   : multiplier for the variance-covariance of
%                                       the active state function
%               |- inactive_var_mult : multiplier for the variance-covariance of
%                                       the inactive state function (default: 0 i.e. a dirac)
%               |- Optim_method      : Choice of optimization routine. By 
%                                      default uses minFunc (Mark Schmidt's implementation). Specify 'fminunc' to use
%                                      Matlab optimization toolbox routine.
%          -clustering       : Structure containing parcelization parameters                        
%               |- MSP_R2_threshold  : threshold (between 0 and 1) to filter 
%                                      data and forward operator in the msp
%                                      calculation.
%               |- MSP_window        : number of samples to use for calculating the MSP. 
%                                      The window should be small enough to satisfy assumption of 
%                                      stationarity within that window.
%               |- clusters_type     : Set to 1 to used C. Grova's clustering
%                                      method
%               |- MSP_scores_threshold : Threshold of the MSP scores (between 0 
%                                         and 1). Set to 0 to keep all the scores.
%               |- neighborhood_order: Nb. of neighbors used to cluterize 
%                                      cortical surface
%          -model          : Parameters for MEM model
%               |-alpha_threshold    : Threshold on the cluster probabilities of
%                                      being active. Clusters whose probas are
%                                      lower than this threshold will be turned off
%               |- alpha             : Option to specify the alpha probabilities 
%                                      manually
%               |- alpha_method      : Selection of the method to initialize the
%                                      probability (alpha) of activation of the 
%                                      clusters  (alpha).
%                       |                   1 = average of the scores of sources in a cluster 
%                       |                   2 = maximum of the scores of sources in a cluster 
%                       |                   3 = median of the scores of sources in a cluster 
%                       |                   4 = equiprobable activity (0.5)
%               |- active_mean_method:  Selection of the method to initialize 
%                                       the mean of the "active state" law (mu).
%                       |                   1 = null hypothesis
%                       |                   2 = Null activity (0)
%          -mandatory          : structure containing data and channel information
%               |- Data              : Data to process (channels x time samples)
%          -automatic         : structure containing outputs
%               |- Units             : Structure for rescaling units
%
% OUTPUTS:
%     - OPTIONS : Keep track of parameters
%     - Results : Structure
%          |- ImageGridAmp  : Source activity (sources x times)
%          |- ImagingKernel : Not computed ([])
%
% ==============================================
% Copyright (C) 2011 - MultiFunkIm & LATIS Team
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO DO LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%    INITIALIZE GLOBAL VARIABLE IF FIRST STUDY
%%%%    CLEAR GLOBAL VARIABLE IF LAST STUDY
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF TO DO LIST %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global MEMglobal

if OPTIONS.optional.verbose
    fprintf('\n\n===== pipeline cMEM\n');
end 

%% Retrieve vertex connectivity - needed for clustering
[OPTIONS, obj.VertConn] = be_vertex_connectivity(HeadModel, OPTIONS);

if isempty(OPTIONS.optional.clustering) && isempty(obj.VertConn) || diff(size(obj.VertConn))
    fprintf('MEM error : no vertex connectivity matrix available.\n');
    return
end

%% ===== Comment ===== %%
OPTIONS.automatic.Comment       =   OPTIONS.optional.Comment;
if strcmp( OPTIONS.automatic.Comment(1:3), 'MEM' )
    OPTIONS.automatic.Comment   =   ['c' OPTIONS.optional.Comment];
end

%% ===== DC offset ===== %% 
% we remove the DC offset the data
[OPTIONS]       = be_remove_dc(OPTIONS);

%% ===== Channels ===== %% 
% we retrieve the channels name and the data
[OPTIONS, obj]  = be_main_channel(HeadModel, obj, OPTIONS);

%% ===== AVG reference ===== %% 
% we average reference the data
[OPTIONS]       = be_avg_reference(OPTIONS);

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

%% ===== Apply temporal data window  ===== %%
% check for a time segment to be localized
[OPTIONS] = be_apply_window( OPTIONS, [] );

%% ===== Double precision to single  ===== %%
% relax the double precision for the msp (leadfield and data)
[OPTIONS] = be_switch_precision( OPTIONS, 'single' );

%% ===== Clusterize cortical surface ===== %%

[OPTIONS, obj] = be_main_clustering(obj, OPTIONS);

%% ===== Single precision to double  ===== %%

[OPTIONS] = be_switch_precision( OPTIONS, 'double' );

%% ===== pre-processing for spatial smoothing (Green mat. square) ===== %%
% matrix W'W from the Henson paper
[OPTIONS, obj.GreenM2] = be_spatial_priorw( OPTIONS, obj.VertConn);


%% ===== Noise estimation ===== %%   

[OPTIONS, obj] = be_main_data_preprocessing(obj, OPTIONS);

%% ===== Fuse modalities ===== %%   

obj = be_fusion_of_modalities( [], obj, OPTIONS);
         	
%% ===== Solve the MEM ===== %%

[obj.ImageGridAmp, OPTIONS] = be_launch_mem(obj, OPTIONS);


%% ===== Inverse temporal data window  ===== %%

[OPTIONS, obj] = be_apply_window( OPTIONS, obj );

% Clean up options
OPTIONS.automatic    = struct(  'entropy_drops', OPTIONS.automatic.entropy_drops, ... 
                                'final_alpha', {OPTIONS.automatic.final_alpha}, ...
                                'Comment', OPTIONS.automatic.Comment, ...
                                'initial_alpha', obj.ALPHA, ...
                                'clusters', obj.CLS);      
                   
% Results (full temporal sequence)                  
Results = struct(...
    'ImageGridAmp',     obj.ImageGridAmp, ...
    'ImagingKernel',    [], ...
    'MEMoptions',       OPTIONS); % ...
     %'MEMdata',          obj);


return




