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
        
        function test_cortical_simulation(testCase)
            
            % Simulate cortical signals -> dHb
            [sHbCortex, beta_map_hb, activation_scout, stim_event_names] = dOD_from_simulated_cortical_activation(testCase.tmp_dir);
            
            % ASSUME: sHbCortex(1) is HbO and sHbCortex(2) is HbR
            % TODO: use common enum for safer identification
            
            % Run GLM
            sGlmResults = cell(1, length(sHbCortex));
            for ihb=1:length(sHbCortex)
                sGlmResults{ihb} = bst_process('CallProcess', 'process_nst_compute_glm', sHbCortex(ihb), [], ...
                                                'stim_events',    strjoin(stim_event_names, ','), ...
                                                'hrf_model',      1, ...  % CANONICAL
                                                'trend',          0, ...
                                                'fitting',        1, ...  % OLS
                                                'save_residuals', 1, ...
                                                'save_betas',     1);
            end
            
            % Check beta estimates, non-regression test at specific voxel
            % where esimtates are the most accurate
            glm_results_hbo = in_bst_matrix(sGlmResults{1}.FileName);
            testCase.assertTrue( abs(beta_map_hb(1, 4708)-glm_results_hbo.beta(1,4708)) < 3e-7);
            glm_results_hbr = in_bst_matrix(sGlmResults{2}.FileName);
            testCase.assertTrue( abs(beta_map_hb(2, 4658)-glm_results_hbr.beta(1,4658)) < 1.5e-7);

            % Check activation detection (p-val thresholding)
            %   
        end
        
    end
    
end

function [hb_cortex, beta_map_hb, activation_scout, stim_event_names] = dOD_from_simulated_cortical_activation(tmp_dir)

%TODO: add second experimental condition

dt = 0.1; %sec
nb_samples = 6000;
time = (0:(nb_samples-1))*dt;
events = db_template('event');
stim_event_names = {'stim'};
events(1).label = stim_event_names{1};
events(1).times = [10 90 120 160 200 260 300 330 360 440 520 550];
events(1).epochs = ones(1, length(events(1).times));
stim_duration = 25; % sec
events(1).times(2,:) = events(1).times(1,:) + stim_duration;

% Insure rounding is consistent
for icond=1:length(events)
    events(icond).samples = round(events(icond).times ./ dt);
    events(icond).times   = events(icond).samples .* dt;
end

hrf_types = process_nst_compute_glm('get_hrf_types');
[X,names] = process_nst_compute_glm('make_design_matrix', time, events, ...
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
activation_scout = bst_create_scout(subject_name, 'cortex', 'activation', seed_vertex_id, 4, 'User scouts');

%% Create functional cortical signals for HbO and HbR

           % Condition1  
delta_hb = [ 0.01        ;... %HbO (mmol.l-1)
            -0.005       ];   %HbR (mmol.l-1)
delta_hb = delta_hb / 1000; % mol.l-1
HB_TYPES = {'HbO', 'HbR'};
surface_file = activation_scout.sSubject.Surface(activation_scout.isurface).FileName;
sCortex = in_tess_bst(surface_file);
nb_vertices = size(sCortex.Vertices,1);
func_map_hb = zeros(length(HB_TYPES), nb_samples, nb_vertices);
beta_map_hb = zeros(length(HB_TYPES), nb_vertices);
for ihb=1:length(HB_TYPES)
    beta_map_hb(ihb, activation_scout.sScout.Vertices) = delta_hb(ihb, 1);
    func_map_hb(ihb, :, :) = X * beta_map_hb(ihb,:) + (rand(nb_samples, nb_vertices) * 0.1 * delta_hb(ihb, 1));
    nst_bst_add_surf_data(squeeze(func_map_hb(ihb,:, :))', time, head_model_holder, ['simu_signal_' HB_TYPES{ihb}], ...
                          [], sStudy, 'simulated', surface_file);
    nst_bst_add_surf_data(beta_map_hb(ihb,:)', 1, head_model_holder, ['simu_beta_' HB_TYPES{ihb}], ...
                          [], sStudy, 'simulated', surface_file);
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

% Add some noise at the channel level
sig_range = max(reordered_dOD(:)) - min(reordered_dOD(:));
reordered_dOD = reordered_dOD + rand(size(reordered_dOD)) * 0.1 * sig_range;

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
hb_cortex = bst_process('CallProcess', ...
                        'process_nst_cortical_projection_mne', sdODInput, []);
                                       
end