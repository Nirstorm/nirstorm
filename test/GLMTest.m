classdef GLMTest < matlab.unittest.TestCase
    
    properties
        tmp_dir
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
            utest_bst_setup();
        end
    end
    
    methods(TestMethodTeardown)
        function tear_down(testCase)
            rmdir(testCase.tmp_dir, 's');
            utest_clean_bst();
        end
    end
    
    methods(Test)
        
        function test_ols_svd(testCase)
            dt = 1; %sec
            nb_samples = 700;
            nb_channels=10;
            nb_active_channels=5;
            
            [X, names,tmp] = design_matrix(0:dt:(nb_samples-1));

            
            [B, active_channels] = create_B_map(nb_channels,nb_active_channels);
            y=X*B + 1e-3*randn(nb_samples,nb_channels);
            
            [B_hat,proj_X] = process_nst_glm_fit('compute_B_svd',y,X);
            testCase.assertTrue( all(all(abs(B - B_hat) < 1e-3)));
        end
        
        function test_AR_ols_fit(testCase)
            dt = 1; %sec
            nb_samples = 700;
            nb_channels=10;
            nb_active_channels=5;
            
            [X, names,tmp] = design_matrix(0:dt:(nb_samples-1));

            
            [B, active_channels] = create_B_map(nb_channels,nb_active_channels);
            y=X*B + 1e-3*randn(nb_samples,nb_channels);
            
            [B_hat, covB, dfe, residuals, mse_residuals] = process_nst_glm_fit('AR1_ols_fit',y, dt, X, 0.1);
            testCase.assertTrue( all(all(abs(B - B_hat) < 1e-3)));  
         end

        function test_channel_simulation_SNR_5(testCase)
            
            method_test=2; % 1 is to test pre-coloring, 2 is to test pre-whitening 
            
            repo_url = nst_get_repository_url();
            data_fns = nst_request_files({{'unittest','nirscomp_data','GraspingRight_Hb_baseline_#1_.nirs'}}, ...
                              1, repo_url);
            nirs_fn = fullfile(testCase.tmp_dir, 'GraspingRight_Hb_baseline_#1_.nirs');
            copyfile(data_fns{1}, nirs_fn);
            
            %% Import data in brainstorm
            [subject_name, sSubject, iSubject] = bst_create_test_subject();
            sDummy = utest_import_nirs_in_bst(nirs_fn, 0);
            
            % Load Noise
            % Todo FIX - resample to 5Hz 
            NewFiles = import_raw_to_db( sDummy.FileName );
            sData = in_bst_data(NewFiles{1});
            
            % extract HbO channels 
            ChannelMat = in_bst_channel(sDummy.ChannelFile);
            ChannelMat.Nirs.Hb={'HbO'};

            iHbO=[];

            i=1;
            for ichan=1:size(ChannelMat.Channel,2)
                if strcmp(ChannelMat.Channel(ichan).Group, 'HbO')
                    channels(i)=ChannelMat.Channel(ichan);
                    iHbO(end+1)=i;
                    i=i+1;
                end    
            end
            ChannelMat.Channel=channels;
            n_chan=size(channels,2);
            ChannelMat.Comment=['NIRS-BRS sensors (' n_chan ')'];
            

            
            
            % create signal 
            [X, names, events] = design_matrix(sData.Time);
            names=names(1:2);
            
            [B, active_channels] = create_B_map(n_chan,round(0.5*n_chan));
            
            signal=X*B;
            noise=sData.F';
            
            SNR=-5; %in dB
            
            y=zeros( size(signal));
            
            for ichann=1:n_chan
                l=get_factor_SNR(signal(:,ichann),noise(:,ichann),SNR);
                disp([' Using l= ' num2str(l) ' for chan ' num2str(ichann) ]);
                y(:,ichann) = signal(:,ichann) + l*noise(:,ichann);
            end

            % Save Data in a new folder 
            iStudy = db_add_condition('test_subject', 'Simulation');
            sStudy = bst_get('Study', iStudy);

            [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
            db_set_channel(iChannelStudy, ChannelMat, 0, 0);

            
            sDataOut = db_template('data');
            sDataOut.F            = y';
            sDataOut.Comment      = 'NIRS dHb';
            sDataOut.ChannelFlag  =  ones(n_chan,1) ;
            sDataOut.Time         = sData.Time;
            sDataOut.DataType     = sData.DataType;
            sDataOut.nAvg         = 1;
            sDataOut.Events       = events;
            sDataOut.DisplayUnits = 'mumol.l-1';
            
            % Generate a new file name in the same folder
            dHbFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_dHb');
            sDataOut.FileName = file_short(dHbFile);
            bst_save(dHbFile, sDataOut, 'v7');
            db_add_data(iStudy, dHbFile, sDataOut);

            sHb = bst_process('GetInputStruct', dHbFile);
            
            
            sHb_filtered = bst_process('CallProcess', 'process_bandpass', sHb, [], ...
                'sensortypes', '', ...
                'highpass',    0.01, ...
                'lowpass',     0.15, ...
                'tranband',    0, ...
                'attenuation', 'strict', ...  % 60dB
                'ver',         '2019', ...  % 2019
                'mirror',      0, ...
                'overwrite',   0);
            
            sGlmResults = bst_process('CallProcess','process_nst_glm_fit', sHb_filtered, [], ...
                                                    'hpf_low_cutoff',         0.01, ...
                                                    'trim_start',             0, ...
                                                    'fitting',                1, ...  % OLS
                                                    'statistical_processing', method_test, ...
                                                    'noise_model',            1, ...  % AR(1)
                                                    'stim_events',            strjoin(names, ','), ...
                                                    'hrf_model',              1, ...  % CANONICAL
                                                    'trend',                  1, ...
                                                    'save_betas',             0, ...
                                                    'save_residuals',         0, ...
                                                    'save_fit',               0);
                
        
             % Check activation detection (p-val thresholding)
            con_sResults_stim1 = bst_process('CallProcess', 'process_nst_glm_contrast', sGlmResults, [], ...
                                                  'Contrast', '[1 -1]');
            stat_sResults_stim1 = bst_process('CallProcess', 'process_nst_glm_contrast_ttest', ...
                                                  con_sResults_stim1, [], ...
                                                  'tail', 'two');
            glm_results_stat = in_bst_matrix(stat_sResults_stim1.FileName);
            thresh_options.pThreshold = 5e-2;
            thresh_options.Correction = 'bonferroni';
            thresh_options.Control = 1;
            
            [glm_results, corr_p] = bst_stat_thresh(glm_results_stat.pmap, thresh_options);
            true_map=zeros(1,n_chan);
            true_map(active_channels)=1;
            
            plot_performance(glm_results_stat,true_map);


            
            [TP,FP,TN,FN]=compute_performance(glm_results, true_map);
            
            %FPR=FP/(FP+TN);
            %testCase.assertLessThan(FPR, 5e-2);
            
            %sensitivity=TP/(TP+FN);
            %testCase.assertGreaterThan(sensitivity,0.8); 
            
        
        end
        
        function test_cortical_simulation(testCase)
            
            % Simulate cortical signals -> dHb
            [sHbCortex, beta_map_hb, activation_scout, stim_event_names] = dOD_from_simulated_cortical_activation(testCase.tmp_dir);
            
            % ASSUME: sHbCortex(1) is HbO and sHbCortex(2) is HbR
            % TODO: use common enum for safer identification
            
            % Run GLM
            sGlmResults = cell(1, length(sHbCortex));
            for ihb=1:length(sHbCortex)
                sGlmResults{ihb} = bst_process('CallProcess', 'process_nst_glm_fit', sHbCortex(ihb), [], ...
                                                'stim_events',    strjoin(stim_event_names, ','), ...
                                                'hrf_model',      1, ...  % CANONICAL
                                                'trend',          1, ...
                                                'fitting',        1, ...  % OLS - precoloring
                                                'save_residuals', 0, ...
                                                'save_betas',     0);
            end
            
            % Check beta estimates, non-regression test at specific voxel
            % where esimtates are the most accurate
            
            % TODO: better handle unit conversion
            %glm_results_hbo = nst_convert_unit(glm_results_hbo, 'mol.l-1', 'mumol.l-1');
            beta_map_hb = beta_map_hb * 1e6; %convert to mumol
            glm_results_hbo = in_bst_matrix(sGlmResults{1}.FileName);
            poi_activ = 4491; % for condition 1
            testCase.assertTrue( abs(beta_map_hb(1, 1, poi_activ)-glm_results_hbo.ImageGridAmp(poi_activ, 1)) < 5);
            poi_activ = 4890; % for condition 2
            testCase.assertTrue( abs(beta_map_hb(1, 2, poi_activ)-glm_results_hbo.ImageGridAmp(poi_activ, 2)) < 1);
            
            glm_results_hbr = in_bst_matrix(sGlmResults{2}.FileName);
            poi_activ = 4539; % for condition 1
            testCase.assertTrue( abs(beta_map_hb(2, 1, poi_activ)-glm_results_hbr.ImageGridAmp(poi_activ, 1)) < 2);
            poi_activ = 4780; % for condition 2
            testCase.assertTrue( abs(beta_map_hb(2, 2, poi_activ)-glm_results_hbr.ImageGridAmp(poi_activ, 2)) < 0.2);

            % Check activation detection (p-val thresholding)
            con_sResults_stim1_hbo = bst_process('CallProcess', 'process_nst_glm_contrast', sGlmResults{1}, [], ...
                                                  'Contrast', '[1 0]');
            stat_sResults_stim1_hbo = bst_process('CallProcess', 'process_nst_glm_contrast_ttest', ...
                                                  con_sResults_stim1_hbo, [], ...
                                                  'tail', 'two');
            con_sResults_stim2_hbo = bst_process('CallProcess', 'process_nst_glm_contrast', sGlmResults{2}, [], ...
                                                  'Contrast', '[0 1]');
            stat_sResults_stim2_hbo = bst_process('CallProcess', 'process_nst_glm_contrast_ttest', ...
                                                  con_sResults_stim2_hbo, [], ...
                                                  'tail', 'two');
                                   
            stat_results_stim1_hbo = in_bst_results(stat_sResults_stim1_hbo.FileName);
            stat_results_stim2_hbo = in_bst_results(stat_sResults_stim2_hbo.FileName);

            thresh_options.pThreshold = 5e-2;
            thresh_options.Correction = 'bonferroni';
            thresh_options.Control = 1;
            
            % Load scalp/cortex geometry and compute depth
            sCortex = in_tess_bst(stat_results_stim1_hbo.SurfaceFile);
            sSubject = bst_get('Subject', sHbCortex.SubjectName);
            sScalp = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
            [scalp_points, depth] = nst_knnsearch(sScalp.Vertices, sCortex.Vertices);
            
            %TODO: fix number of tests according to cortical coverage
            [pmask_stim1_hbo, corr_p] = bst_stat_thresh(stat_results_stim1_hbo.pmap, thresh_options);
            [pmask_stim2_hbo, corr_p] = bst_stat_thresh(stat_results_stim2_hbo.pmap, thresh_options);

            activation_mask_stim1_hbo = zeros(size(pmask_stim1_hbo));
            activation_mask_stim1_hbo(activation_scout(1).sScout.Vertices) = 1;

            activation_mask_stim2_hbo = zeros(size(pmask_stim2_hbo));
            activation_mask_stim2_hbo(activation_scout(2).sScout.Vertices) = 1;
            
            % Test that all truely activated vertices are detected 
            % limit to vertices not too deep (<3cm)
            %TODO: check activation tests again
            max_depth = 0.03; % meter
            testCase.assertTrue(all(pmask_stim1_hbo(activation_mask_stim1_hbo & depth<max_depth)==1));
            testCase.assertTrue(all(pmask_stim2_hbo(activation_mask_stim2_hbo & depth<max_depth)==1));

            % Non regression test on the number of false positive
            % TODO: improve / fix
            if sum(pmask_stim1_hbo(~activation_mask_stim1_hbo)) > 300
                warning('GLM: number of false positives too large for stim1');
            end
            if sum(pmask_stim2_hbo(~activation_mask_stim2_hbo)) > 1
                warning('GLM: number of false positives too large for stim2');
            end
            
%             testCase.assertLessThan(sum(pmask_stim1_hbo(~activation_mask_stim1_hbo)), 300);
%             testCase.assertLessThan(sum(pmask_stim2_hbo(~activation_mask_stim2_hbo)), 1);

            % TODO: test bilateral            
        end
        
    end
    
end

function converted = nst_convert_unit(values, unit_in, unit_out)
%TODO: test & doc
assert(units_homogeneous(unit_in, unit_out));
converted = values/get_unit_factor(unit_in) * get_unit_factor(unit_out);
end

function [hb_cortex, beta_map_hb, activation_scout, stim_event_names] = dOD_from_simulated_cortical_activation(tmp_dir)

% Returned beta_map_hb is in mol.l-1

dt = 0.1; %sec
nb_samples = 6000;
time = (0:(nb_samples-1))*dt;
events = db_template('event');
stim_event_names = {'stim1', 'stim2'};
events(1).label = stim_event_names{1};
events(1).times = [10 160 200 260 300 330 360 400 440 460 490 520 550];
events(1).epochs = ones(1, length(events(1).times));

events(2).label = stim_event_names{2};
events(2).times = [40 60 90 120 140 150 170 210 270 300 380];
events(2).epochs = ones(1, length(events(2).times));

stim_duration = 10; % sec
events(1).times(2,:) = events(1).times(1,:) + stim_duration;
events(2).times(2,:) = events(2).times(1,:) + stim_duration;

% Insure rounding is consistent
for icond=1:length(events)
    events(icond).samples = round(events(icond).times ./ dt);
    events(icond).times   = events(icond).samples .* dt;
end

hrf_types = process_nst_glm_fit('get_hrf_types');
[X,names] = process_nst_glm_fit('make_design_matrix', time, events, ...
                                hrf_types.CANONICAL, 25, 0);
assert(all(strcmp(names, {events.label})));


%% Retrieve data
repo_url = nst_get_repository_url();
data_fns = nst_request_files({{'unittest','lesca_data','dummy_frontal_16x16.nirs'}, ...
                              {'unittest','lesca_data','optodes_frontal_16x16_Colin27_4NIRS.txt'}, ...
                              {'unittest','lesca_data','headmodel_optodes_frontal_16x16_Colin27_4NIRS.mat'}}, ...
                              1, repo_url);
nirs_fn = fullfile(tmp_dir, 'dummy_frontal_16x16.nirs');
copyfile(data_fns{1}, nirs_fn);
copyfile(data_fns{2}, fullfile(tmp_dir, 'optodes.txt'));

%% Import data in brainstorm
[subject_name, sSubject, iSubject] = bst_create_test_subject();
% Use lowres mid as default surface: 
db_surface_default(iSubject, 'Cortex', find(strcmp({sSubject.Surface.Comment}, 'mid_lowres'))); 
sDummy = utest_import_nirs_in_bst(nirs_fn, 0);

%% Inject headmodel
[condition_dir, bfn, ext] = fileparts(file_fullpath(sDummy.FileName));
headmodel_fn = fullfile(condition_dir, 'headmodel_nirs_mcx_fluence.mat');
copyfile(data_fns{3},  headmodel_fn);

sStudy = bst_get('Study', sDummy.iStudy);
sStudy.iHeadModel = 1;
head_model_holder = db_template('Headmodel');
head_model_holder.FileName = file_short(headmodel_fn);
head_model_holder.Comment = 'NIRS head model (all pairs)';
head_model_holder.HeadModelType = 'surface';
sStudy.HeadModel(1) = head_model_holder;
bst_set('Study', sDummy.iStudy, sStudy);
db_save();

%% Create scout for activating region
seed_vertex_id = 4560; %scout seed 
activation_scout(1) = bst_create_scout(subject_name, 'cortex', 'activation_stim1', seed_vertex_id, 4, 'User scouts');

seed_vertex_id = 4823; %scout seed 
activation_scout(2) = bst_create_scout(subject_name, 'cortex', 'activation_stim2', seed_vertex_id, 3, 'User scouts');

%% Create functional cortical signals for HbO and HbR

           % stim1   stim2  
delta_hb = [ 0.01    0.005     ;... %HbO (mmol.l-1)
            -0.005  -0.002    ];   %HbR (mmol.l-1)
delta_hb = delta_hb / 1000; % mol.l-1
HB_TYPES = {'HbO', 'HbR'};
surface_file = activation_scout(1).sSubject.Surface(activation_scout(1).isurface).FileName;
sCortex = in_tess_bst(surface_file);
nb_vertices = size(sCortex.Vertices,1);
nb_conditions = length(events);
func_map_hb = zeros(length(HB_TYPES), nb_samples, nb_vertices);
beta_map_hb = zeros(length(HB_TYPES), nb_conditions, nb_vertices);
for ihb=1:length(HB_TYPES)
    for icond=1:nb_conditions
        beta_map_hb(ihb, icond, activation_scout(icond).sScout.Vertices) = delta_hb(ihb, icond);
        nst_bst_add_surf_data(squeeze(beta_map_hb(ihb, icond, :)), 1, head_model_holder, ['simu_beta_' HB_TYPES{ihb} '_' events(icond).label], ...
                              ['Simulated beta ' HB_TYPES{ihb} ' ' events(icond).label], [], sStudy, 'simulated', surface_file);
    end
    func_map_hb(ihb, :, :) = X * squeeze(beta_map_hb(ihb,:,:) ) + rand(nb_samples, nb_vertices) * 0.1 * mean(delta_hb(ihb,:));
    nst_bst_add_surf_data(squeeze(func_map_hb(ihb,:, :))', time, head_model_holder, ['simu_sig_' HB_TYPES{ihb}], ...
                          ['Simulated signals ' HB_TYPES{ihb}], [], sStudy, 'simulated', surface_file);
end
%% Project in topographic space using head model
head_model = in_bst_headmodel(headmodel_fn);
channel_def = in_bst_channel(sDummy.ChannelFile);
montage_info = nst_montage_info_from_bst_channels(channel_def.Channel);
sensitivity_surf = process_nst_import_head_model('get_sensitivity_from_chans', head_model, montage_info.pair_names);

A_w1 = squeeze(sensitivity_surf(:,1,:));
A_w2 = squeeze(sensitivity_surf(:,2,:));

% nb_wavelengths x [HbO, HbR]
ext_coeffs = process_nst_mbll('get_hb_extinctions', channel_def.Nirs.Wavelengths);

g = [A_w1*ext_coeffs(1,1) A_w1*ext_coeffs(1,2);...
     A_w2*ext_coeffs(2,1) A_w2*ext_coeffs(2,2)];
func_map_hb_stacked = [squeeze(func_map_hb(1,:,:))' ; squeeze(func_map_hb(2,:,:))'];

func_map_dOD_stacked = g * func_map_hb_stacked;

% Align functional data to current channel definition
reordered_dOD = zeros(size(func_map_dOD_stacked));
nb_pairs = size(montage_info.pair_ichans, 1);
for ipair=1:nb_pairs
    for imeasure=1:size(montage_info.pair_ichans, 2)
        ichan = montage_info.pair_ichans(ipair, imeasure);
        reordered_dOD(ichan, :) = func_map_dOD_stacked(ipair + nb_pairs*(imeasure-1), :);
    end
end

% Add some white noise at the channel level
sig_range = max(reordered_dOD(:)) - min(reordered_dOD(:));
reordered_dOD = reordered_dOD + randn(size(reordered_dOD)) * 0.1 * sig_range;

%% Save dOD as new data
sDummyData = in_bst_data(sDummy.FileName);
sDataOut = db_template('data');
sDataOut.F            = reordered_dOD;
sDataOut.Comment      = 'NIRS dOD simulated';
sDataOut.ChannelFlag  = sDummyData.ChannelFlag;
sDataOut.Time         = time;
sDataOut.DataType     = 'recordings';
sDataOut.nAvg         = 1;
sDataOut.Events       = events;
sDataOut.DisplayUnits = 'delta OD';

% Generate a new file name in the same folder
dODFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_dod');
sDataOut.FileName = file_short(dODFile);
bst_save(dODFile, sDataOut, 'v7');
% Register in database
db_add_data(sDummy.iStudy, dODFile, sDataOut);

%TODO: check channels that are on top of activation area but no evoked signal

%% Project back on the cortex

sdODInput = bst_process('GetInputStruct', dODFile);
proj_methods = process_nst_cortical_projection('methods');
hb_cortex = bst_process('CallProcess', ...
                        'process_nst_cortical_projection', sdODInput, [], ...
                        'method', proj_methods.MNE);                  
end


function [X, names, events] = design_matrix(time)

% Returned the design matrix
% generate events with the repetition of the following paterns 
% Stim 1 [ 20 - 30s]  - Stim 2 [20-30s] 


events = db_template('event');
stim_event_names = {'stim1', 'stim2'};

events(1).label = stim_event_names{1};
events(2).label = stim_event_names{2};

t=0;
i_event=1;

while t + 60 < time(end)
    
    duration_stim1=randi([20,30]);
    duration_stim2=randi([20,30]);
    
    events(1).times(1,i_event)=t;
    events(1).times(2,i_event)=t+duration_stim1;
    
    events(2).times(1,i_event)=t+duration_stim1;
    events(2).times(2,i_event)=t+duration_stim1+duration_stim2;
   
    t=t+duration_stim1+duration_stim2;
    i_event=i_event+1;
    
end

events(1).epochs = ones(1, length(events(1).times));
events(2).epochs = ones(1, length(events(2).times));

hrf_types = process_nst_glm_fit('get_hrf_types');
[X,names] = process_nst_glm_fit('make_design_matrix', time, events, ...
                                hrf_types.CANONICAL, 25, 1);

                            
assert(all(strcmp(names(1:2), {events.label})));

end

function [B, active_channels]= create_B_map(nb_channels,nb_active_channels)
    % Assume 2 conditions Stim1=Task,Stim2=Rest  
    B=zeros(3,nb_channels);
    perm=randperm(nb_channels);
    
    active_channels=sort(perm(1:nb_active_channels));
    unactive_channels=sort(perm(nb_active_channels+1:end));
    
    
    % Active channels : B1 = 0.01; B2~ 0.005, B3=N(0,0.0001)
    B(1,active_channels)=0.01+0.0001*randn(); % Stim 1
    B(2,active_channels)=0.005+0.0001*randn(); % Stim 2 
    B(3,active_channels)=0.0001*randn(); % Constant 

    % Active channels : B1 =  N(0,0.001); B2~ N(0,0.001), B3=0

    B(1,unactive_channels)=0.001*randn();% Stim 1
    B(2,unactive_channels)=0.001*randn(); % Stim 2 
    B(3,unactive_channels)=0.0001*randn(); % Constant 
    
    B=B/100;
end


function l=get_factor_SNR(signal,noise,SNR)
    % compute l so that y = signal + l*noise avec the correct SNR value
    % SNR = 10 log(  <signal> / l^2 < noise > ) with <.> the signal power
    
    mss_signal=norm(signal,2);
    mss_noise=norm(noise,2);
    
    l2= (mss_noise/mss_signal) * 10^(SNR/10);
    l=sqrt(l2);
end

function [TP,FP,TN,FN]=compute_performance(glm_results, true_results)
    FP= sum(glm_results( true_results == 0  ));
    TP= sum(glm_results( true_results >  0 ));
    
    TN= sum(~glm_results(true_results == 0));
    FN= sum(~glm_results(true_results  >  0 ));
end

function plot_performance(glm_results_stat,true_map)

    pval= logspace(-7,0,1000);
    
    N=length(pval);
    
    FPR=zeros(1,N);
    sensitivity=zeros(1,N);

    for i=1:N
        thresh_options.pThreshold = pval(i);
        thresh_options.Correction = 'bonferroni';
        thresh_options.Control = 1;

        [glm_results, corr_p] = bst_stat_thresh(glm_results_stat.pmap, thresh_options);
        [TP,FP,TN,FN]=compute_performance(glm_results, true_map);

        FPR(i)=FP/(FN+TP);
        sensitivity(i)=TP/(TP+FN);
        
    end
    
    figure(1)
    plot(FPR,sensitivity);
    xlabel('False Positive Rate')
    ylabel('Sensitivity')
    
    figure(2)
    plot(pval,FPR);
    xlabel('p-value (Bonferroni correction)')
    ylabel('False Positive Rate')
    
    figure(3)
    plot(pval,sensitivity);
    xlabel('p-value (Bonferroni correction)')
    ylabel('Sensitivity')

    save( ['/NAS/home/edelaire/Documents/output_analysis/simulation/GLM_whit_filterd_FPR.mat'], 'FPR');
    save( ['/NAS/home/edelaire/Documents/output_analysis/simulation/GLM_whit_filterd_sens.mat'], 'sensitivity');

    
end