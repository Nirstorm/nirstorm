function [Results, OPTIONS] = be_rwmem_solver(HeadModel, OPTIONS, Results)
% RWMEMSOLVER: Maximum Entropy on the Mean solution.
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


% pipeline starts here:
fprintf('\n\n===== pipeline rwMEM\n');
        
        
%% Useful variables
obj.ImageGridAmp = []; 
if ~isfield(HeadModel, 'vertex_connectivity')
    [OPTIONS, obj.VertConn] = be_vertex_connectivity(HeadModel, OPTIONS);
else
    obj.VertConn = HeadModel.vertex_connectivity;
end

if isempty(OPTIONS.optional.clustering) && ( isempty(obj.VertConn) || diff(size(obj.VertConn)) ) % CETTE PARTIE N'EST PAS CORRECTE: || size(obj.VertConn,1)~=size(HeadModel.Gain,2)
    fprintf('MEM error : no vertex con.nectivity matrix available.\n');
    return
end


% Get filtered baseline
DISP    =   OPTIONS.optional.display;
if ~OPTIONS.automatic.stand_alone
    TFf     =   load( be_fullfile(OPTIONS.automatic.iProtocol.STUDIES, OPTIONS.automatic.DataInfo.TFfile), 'Options' );
else
    TFf     =   load( OPTIONS.automatic.DataInfo.TFfile, 'Options' );
end
BSL     =   be_bpfilter( OPTIONS.automatic.DataInfo.Baseline, OPTIONS.automatic.sampling_rate, TFf.Options.RidgesFRange );


% Real part
fprintf('rwMEM on analytic signal...\treal part...\t')
OPTIONS.mandatory.pipeline  =   'wMEM';
OPTIONS.optional.Baseline   =   BSL;
OPTIONS.optional.display    =   0;
[RealR]                     =   be_wmem_solver(HeadModel, OPTIONS, Results);
if DISP
    [Hr.hfig Hr.hfigtab]    =   be_display_time_scale_boxes(RealR.MEMdata,RealR.MEMoptions);
    set(Hr.hfig, 'Name', 'Real analytic signal');
    be_display_entropy_drops(Hr,RealR.MEMoptions);    
end
fprintf('Done.\n');


% Get imaginary signal 
TMbin   =   be_closest(OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime);
OPTIONS.mandatory.Data      =   OPTIONS.automatic.DataInfo.Fi( :, TMbin(1):TMbin(end) );   

% Imaginary part
fprintf('rwMEM on analytic signal...\timaginary part...\t')
[ImagR]                     =   be_wmem_solver(HeadModel, OPTIONS, Results);
if DISP    
    [Hi.hfig Hi.hfigtab]    =   be_display_time_scale_boxes(ImagR.MEMdata,ImagR.MEMoptions);
    set(Hi.hfig, 'Name', 'Imaginary analytic signal');
    be_display_entropy_drops(Hi,ImagR.MEMoptions);
end
fprintf('Done.\n');


% Organize figure;
if DISP 
    set(Hr.hfig, 'Units', 'Normalized', 'Position', [.1 .3 .4 .5]);
    set(Hi.hfig, 'Units', 'Normalized', 'Position', [.5 .3 .4 .5]);
end


% Results (full temporal sequence)
Results = struct(...
    'ImageGridAmp',     RealR.ImageGridAmp+1i*ImagR.ImageGridAmp, ...
    'ImagingKernel',    [], ...
    'MEMoptions',       fuse(RealR.MEMoptions, ImagR.MEMoptions,1), ...
    'MEMdata',          fuse(RealR.MEMdata, ImagR.MEMdata,2) );

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- HELPER FUNCTIONS --------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [OUT] = fuse(Real, Imag, what)

    switch what
        
        case 1      % OPTIONS
            
            % Output
            OUT                             =   Real;
            
            % Automatic
            OUT.automatic                   =   Real.automatic;
            OUT.automatic.Modality(2)       =   Imag.automatic.Modality;
            OUT.automatic.final_alpha       =   {OUT.automatic.final_alpha Imag.automatic.final_alpha};
            OUT.automatic.entropy_drops     =   [{OUT.automatic.entropy_drops} {Imag.automatic.entropy_drops}];
            OUT.automatic.selected_samples  =   {OUT.automatic.selected_samples Imag.automatic.selected_samples};
            OUT.automatic.selected_values   =   {OUT.automatic.selected_values Imag.automatic.selected_values};
            OUT.automatic.Mod_in_boxes      =   {OUT.automatic.Mod_in_boxes Imag.automatic.Mod_in_boxes};
            OUT.automatic.scales            =   {OUT.automatic.scales Imag.automatic.scales};
            OUT.automatic.final_sigma       =   {OUT.automatic.final_sigma Imag.automatic.final_sigma};
            
            % Optional
            OUT.optional.iData              =   Imag.mandatory.Data;
            SEPcomm                         =   strfind( OUT.automatic.DataInfo.Comment, '|' );
            OUT.automatic.Comment           =   ['r' OUT.optional.Comment OUT.automatic.DataInfo.Comment(SEPcomm-1:end)];
            
        case 2      % OBJ
            
            % Output
            OUT             =   Real;
            OUT             =   rmfield( OUT, 'ImageGridAmp' );
            OUT.data        =   {OUT.data Imag.data};
            OUT.SCR         =   {OUT.SCR Imag.SCR};
            OUT.CLS         =   {OUT.CLS Imag.CLS};
            OUT.ALPHA       =   {OUT.ALPHA Imag.ALPHA};
            OUT.noise_var   =   {OUT.noise_var Imag.noise_var};
            OUT.baseline    =   {OUT.baseline Imag.baseline};
            
    end
            

return


