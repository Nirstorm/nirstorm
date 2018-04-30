function [Results, OPTIONS] = be_rmem_solver(HeadModel, OPTIONS, Results)
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


% Display progress
if OPTIONS.optional.waitbar
    obj.hmem = waitbar(0, {'Step 1/2 : Preparing MEM routine ...'});
end


%% Retrieve vertex connectivity - needed for clustering
if ~OPTIONS.automatic.stand_alone
    [OPTIONS, obj.VertConn] = be_vertex_connectivity(HeadModel, OPTIONS);
elseif isempty(OPTIONS.optional.clustering.clusters) % we consider a stand-alone version (outside the BrainStorm env.) 
    obj.VertConn = HeadModel.vertex_connectivity;
end
obj.ImageGridAmp = []; 

if isempty(OPTIONS.optional.clustering) && isempty(obj.VertConn) || diff(size(obj.VertConn)) 
    fprintf('MEM error : no vertex connectivity matrix available.\n');
    return
end

% if OPTIONS.optional.waitbar
%     % Display progress
%     waitbar(0.1, obj.hmem, {'Step 1/2 : Preparing MEM routine ... head files loaded'});
% end

%% ===== Comment ===== %%
OPTIONS.automatic.Comment       =   OPTIONS.optional.Comment;
if strcmp( OPTIONS.automatic.Comment(1:3), 'MEM' )
    OPTIONS.automatic.Comment   =   ['r' OPTIONS.optional.Comment];
end

%% ===== Channels ===== %% 

[OPTIONS, obj] = be_main_channel(HeadModel, obj, OPTIONS);

%% ===== Sources ===== %% 

[OPTIONS, obj] = be_main_sources(obj, OPTIONS);

%% ===== Pre-process the leadfield(s) ===== %% 

[OPTIONS, obj] = be_main_leadfields(obj, OPTIONS);


%% ===== Compute ridge-filtered signals ===== %%
    
[OPTIONS, obj] = get_ridge_data(OPTIONS, obj, []);


%% ===== Process NEW ridge-filtered signals ==== %%

nNEWdata    =   numel(OPTIONS.automatic.rMEMfiles);
Results     =   struct( 'ImageGridAmp',     [], ...
                        'ImagingKernel',    [], ...
                        'MEMoptions',       OPTIONS, ...
                        'MEMdata',          obj);
    
for ii = 1 : nNEWdata
    
    fprintf('\nProcessing ridge #%i ...\n', ii)
    
    %% ===== Load Data ===== %%
    
    [OPTIONS2, obj] = get_ridge_data(OPTIONS, obj, ii);
    
    %% ===== Normalization ==== %% 
    % we absorb units (pT, nA) in the data, leadfields; we normalize the data
    % and the leadfields
    [OPTIONS2] = be_normalize_and_units(OPTIONS2);

    %% ===== Noise estimation ===== %%   

    [OPTIONS2, obj] = be_main_data_preprocessing(obj, OPTIONS2);

    %% ===== Clusterize cortical surface ===== %%
    
    [OPTIONS2, obj] = be_main_clustering(obj, OPTIONS2);
    
    %% ===== pre-processing for spatial smoothing (Green mat. square) ===== %%
    % matrix W'W from the Henson paper
    [OPTIONS2, obj.GreenM2] = be_spatial_priorw( OPTIONS2, obj.VertConn);
    
    if OPTIONS.optional.waitbar
        % Display progress
        waitbar(0.95, obj.hmem, {'Step 1/2 : Preparing MEM routine ... spatial smoothing operator computed'});
    end
    
    %% ===== double precision ===== %%
    [OPTIONS2]       = be_switch_precision( OPTIONS2, 'double' );
        
    %% ===== Solve the MEM ===== %%
    
    [obj, OPTIONS2] = be_main_rmem(obj, OPTIONS2);
    
    % new filename
    [rep, base, ext] = be_fileparts(OPTIONS2.optional.ResultFile);
    OPTIONS2.optional.ResultFile = be_fullfile(rep, [base, '_', num2str(ii), ext]);
    
    if OPTIONS2.automatic.stand_alone & ~OPTIONS2.automatic.process
        Results(ii).ImageGridAmp    =   obj.ImageGridAmp;
        Results(ii).ImagingKernel   =   [];
        Results(ii).MEMoptions      =   OPTIONS2;
        Results(ii).MEMdata         =   rmfield(obj, 'ImageGridAmp');
    elseif ii<nNEWdata 
        [OPTIONS2] = save_ridges_results(obj, OPTIONS2, HeadModel, ii);
    else
        [OPTIONS2, Results] = save_ridges_results(obj, OPTIONS2, HeadModel, ii);
    end
end

if OPTIONS.optional.waitbar
    pause(.5); close(obj.hmem)
end

if ~OPTIONS.automatic.stand_alone | OPTIONS.automatic.process                                 
	
 	% If added to a 'default_study' node: need to update results links12 30
	OPTIONS.automatic           =   OPTIONS2.automatic;
    OPTIONS.mandatory.DataTime  =   Results.Time;
 	OPTIONS.optional.DataFile   =   OPTIONS.automatic.rMEMfiles{end};
    
    fprintf('\n***\tWARNING\t***\n\tTo see the ridges and localizations, you need to reload the condition\n\tRight-click on condition -> Reload (at the bottom of the menu)\n\n');
end

return





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- HELPER FUNCTIONS --------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [OPTIONS, Results] = save_ridges_results(obj, OPTIONS, HM, ii)
% ===== SAVE RESULTS FILE =====
    
    % Add a few fields (for records)
    iP                      = bst_get('ProtocolInfo');
    Results.ImageGridAmp    = obj.ImageGridAmp;
    Results.ImagingKernel   = [];    
    Results.nComponents     = round( size(obj.ImageGridAmp,1) / obj.nb_sources );
    Results.Comment         = OPTIONS.automatic.Comment;
    Results.Function        = OPTIONS.InverseMethod;
    Results.Pipeline        = OPTIONS.mandatory.pipeline;
    Results.Time            = load( be_fullfile(OPTIONS.automatic.iProtocol.STUDIES, OPTIONS.automatic.rMEMfiles{ii}), 'Time' );
    Results.Time            = Results.Time.Time;
    Results.ChannelFlag     = OPTIONS.optional.ChannelFlag;
    Results.GoodChannel     = OPTIONS.automatic.GoodChannel;
    Results.HeadModelFile   = file_win2unix(strrep(OPTIONS.optional.HeadModelFile, iP.STUDIES, ''));
    Results.HeadModelType   = HM.HeadModelType;
    Results.SurfaceFile     = file_win2unix(strrep(HM.SurfaceFile, iP.SUBJECTS, ''));
    Results.DataFile        = file_win2unix(strrep(OPTIONS.automatic.rMEMfiles{ii}, iP.STUDIES, ''));
    Results                 = bst_history('add', Results, 'compute', ['Source estimation: ' OPTIONS.InverseMethod ', Pipeline : ' OPTIONS.mandatory.pipeline]);
    %obj = rmfield(obj, 'ImageGridAmp');
    %Results.MEMdata         = obj;
    Results.MEMoptions      = OPTIONS;
    Results.MEMoptions      = setfield(Results.MEMoptions, 'automatic', ...
                                struct(  'entropy_drops', OPTIONS.automatic.entropy_drops, ... 
                                'final_alpha', {OPTIONS.automatic.final_alpha}, ...
                                'Comment', OPTIONS.automatic.Comment, ...
                                'initial_alpha', obj.ALPHA, ...
                                'clusters', obj.CLS) );    
    
    if strcmpi(Results.HeadModelType, 'volume')
        Results.GridLoc     = HeadModel.GridLoc;
    else
        Results.GridLoc     = [];
    end
    
    % Save the relevant method options
    strOptions = fieldnames(OPTIONS);
    for optField = strOptions
        Results.Options.(optField{1}) = OPTIONS.(optField{1});
    end
    
    if ii==numel(OPTIONS.automatic.rMEMfiles)
        return
    end
    
    % Save Results structure
    save( be_fullfile(iP.STUDIES, OPTIONS.optional.ResultFile), '-struct', 'Results');

    % Add a new Result entry to the Brainstorm database
    iStudy                  = be_get_id(OPTIONS);
    sStudy                  = bst_get('Study', iStudy );
    newResult               = db_template('Results');
    newResult(1).Comment    = Results.Comment;
    newResult.FileName      = file_win2unix(strrep(OPTIONS.optional.ResultFile, iP.STUDIES, ''));
    newResult.DataFile      = Results.DataFile;
    newResult.isLink        = 0;
    newResult.HeadModelType = Results.HeadModelType;
    % Get relative Data filename
    newResult.DataFile      = strrep(newResult.DataFile, iP.STUDIES, '');
    % Update Brainstorm DataBase
    iResult                 = length(sStudy.Result) + 1;
    sStudy.Result(iResult)  = newResult;
    bst_set('Study', iStudy, sStudy);
    db_save();    
              
return

function [OPTIONS, obj] = get_ridge_data(OPTIONS, obj, ii)

    iP  =   OPTIONS.automatic.iProtocol;    
    
    % Retrieve ridges infos
    if isempty(ii) && ~isempty(OPTIONS.automatic.rMEMfiles)
        [obj]           = be_fusion_of_modalities([], obj, OPTIONS);
        
    elseif isempty(ii) && isempty(OPTIONS.automatic.rMEMfiles)
        [obj]           = be_fusion_of_modalities([], obj, OPTIONS);
        
        % single precision
        [OPTIONS]       = be_switch_precision( OPTIONS, 'single' );
        
        [OPTIONS, rF]   = be_ridgefilter(OPTIONS);
        
        OPTIONS.automatic.rMEMfiles = rF;
        obj.CLS         = [];
        obj.ALPHA       = [];
        
    elseif ~isempty(ii)
        
        DT              = OPTIONS.automatic.rMEMfiles{ii};
        
        if ischar(DT)
            DT              =   load( be_fullfile(iP.STUDIES, OPTIONS.automatic.rMEMfiles{ii}) );
            FRQs            =   load( be_fullfile(iP.STUDIES, DT.TFfile), 'Freqs', 'Options', 'TFpath', 'TFfile');
            FRQr            =   FRQs.Options.RidgesFRange;
            mspDATA         =   load( be_fullfile( iP.STUDIES, DT.History{end} ), 'F', 'Time');
            mspDATA.FRQs    =   FRQs.Freqs;        
        
        elseif isstruct(DT)
            FRQs.Freqs      =   DT.TFfile.Freqs;
            FRQr            =   DT.TFfile.Options.RidgesFRange;
            mspDATA.F       =   zeros( numel(OPTIONS.mandatory.ChannelTypes), numel(OPTIONS.mandatory.DataTime) );
            mspDATA.Time    =   OPTIONS.mandatory.DataTime;
            mspDATA.FRQs    =   DT.TFfile.Freqs;
            mspDATA.F(OPTIONS.automatic.GoodChannel,: ) = obj.data;
        
        end
        
        % Fill optional fields : DATA
        iDc     = OPTIONS.automatic.GoodChannel;
        STd     = be_closest( OPTIONS.optional.TimeSegment(1), DT.Time );
        NDd     = be_closest( OPTIONS.optional.TimeSegment(end), DT.Time );
        OPTIONS.mandatory.Data              =   DT.F( iDc, STd:NDd );
        OPTIONS.optional.iData              =   DT.Fi( iDc, STd:NDd );  
        
        % Fill optional field : BASELINE
        OPTIONS.optional.Baseline           =   DT.Baseline(iDc, :);
        noise_covariance                    =   [];
        switch OPTIONS.solver.NoiseCov_method
            
            case {4,5}
                % Nothing to do
                % OPTIONS.optional.Baseline(iDc,:)    =   be_bpfilter( OPTIONS.optional.Baseline(iDc,:), OPTIONS.automatic.sampling_rate, FRQr );
                
            case {1, 2, 3}
                FRQidx                      =   mod(DT.TFpath, numel(FRQs.Freqs));
                FRQidx(FRQidx==0)          	=   numel(FRQs.Freqs);
                O                        	=   OPTIONS;
                O.ridges.frequency_range   	=   FRQs.Freqs( ceil(median(FRQidx)) ) + [-3 3];
                noise_covariance          	=   be_randomize_ridge_baseline(DT.Baseline(iDc,:), O);          
        end
        
        % Fill optional fields : TIME DEFINITIONS
        OPTIONS.mandatory.DataTime          =   DT.Time(STd:NDd);
        OPTIONS.ridges.fdrRdgMod            =   DT.rdgMOD;
        OPTIONS.automatic.Comment           =   [OPTIONS.optional.Comment ' ridgefiltered signal # ' num2str(ii)];
        OPTIONS.wavelet.freqs_analyzed      =   FRQs.Freqs;
        OPTIONS.automatic.selected_samples  =   DT.TFpath;
        
        % Filter MSP DATA
        STm                 =   be_closest( OPTIONS.optional.TimeSegment(1), mspDATA.Time );
        NDm                 =   be_closest( OPTIONS.optional.TimeSegment(end), mspDATA.Time );
        %temp_sig            =   mspDATA.F(iDc, :);
        mspDATA.F           =   be_bpfilter( mspDATA.F(iDc,:), OPTIONS.automatic.sampling_rate, FRQr );
        if strcmp(OPTIONS.clustering.clusters_type, 'blockwise')
            mspDATA.F       =   mspDATA.F(:,STm:NDm);
            %temp_sig        =   temp_sig(:,STm:NDm);
        end
        OPTIONS.temporary.mspDATA =   mspDATA;
        
        % ADJUST BASELINE AMPLITUDE
%         sigA                =   std( temp_sig(:) );        
%         tempDT              =   OPTIONS.mandatory.Data(OPTIONS.automatic.GoodChannel,~any(isnan(OPTIONS.mandatory.Data(OPTIONS.automatic.GoodChannel,:)))); 
%         tempDTi             =   OPTIONS.optional.iData(OPTIONS.automatic.GoodChannel,~any(isnan(OPTIONS.optional.iData(OPTIONS.automatic.GoodChannel,:))));
%         rdgR                =   std( tempDT(:) );
%         rdgI                =   std( tempDTi(:) );
%         RATIO               =   mean([rdgR rdgI]) / sigA;
%         OPTIONS.mandatory.Data  =   OPTIONS.mandatory.Data / RATIO;
%         OPTIONS.optional.iData  =   OPTIONS.optional.iData / RATIO;        
        
        % CREATE NEW MODALITY FIELD IN OPTIONS.AUTOMATIC
        [OPTIONS, obj]  =   be_main_channel(0, obj, OPTIONS);
        OPTIONS         =   rmfield( OPTIONS, 'temporary' );
        OPTIONS.automatic.Modality.covariance   =   [];
        OPTIONS.solver.NoiseCov                 =   noise_covariance;
        OPTIONS.solver.NoiseCov_recompute       =   0;
        
    end
    
return

