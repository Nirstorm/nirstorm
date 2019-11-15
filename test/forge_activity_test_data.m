function forge_activity_test_data(dest_dir)
% Generate test data for activity detection processes.
% The generated file is a bst subject that is intended to be uploaded
% on the nirstorm data repository. Then this data is retrieved by
% nst_request_data({{'unittest', 'activity_test_subject', 'activity_test_subject.zip'},...
%                   {'unittest', 'activity_test_subject', 'activity_info.mat'}});
% So this function should not be used directly within unit tests.
% See load_activity_test_subject
%
% The anatomy is Colin27_4NIRS.
%
% The generated data contains three separate regions of interest 
% with no signal overlap.
% Distinct activation patterns are:
%  - frontal activates for stim1 but not stim2
%  - occipital activates for stim2 but not stim1
%  - top does not activate (only noise)
%
% Channel-space and cortical signals are available.
% However there were not built using a proper head model, so there is no
% scaling consistency between them. 
% The intended usage is to analyse them separately.

if nargin < 1
    dest_dir = fullfile(nst_get_local_user_dir(), 'unittest', ...   
                        'activity_test_subject');
    if ~exist(dest_dir, 'dir')
        mkdir(dest_dir);
    end
end
assert(exist(dest_dir, 'dir')~=0)

utest_bst_setup();

bfn = 'test_data.nirs';
nirs_data_fn = fullfile(dest_dir, bfn);

% Time definition
dt = 0.1; %sec
nb_samples = 2400;
time = (0:(nb_samples-1))*dt;

%% Paradigm
events = db_template('event');
stim_event_names = {'stim1', 'stim2'};
events(1).label = stim_event_names{1};
events(1).times = [10 40 66 120 150 182];
events(1).epochs = ones(1, length(events(1).times));

events(2).label = stim_event_names{2};
events(2).times = [20 51 90 110 125 160 190 210];
events(2).epochs = ones(1, length(events(2).times));

stim_duration = 10; % sec
events(1).times(2,:) = events(1).times(1,:) + stim_duration;
events(2).times(2,:) = events(2).times(1,:) + stim_duration;

% Insure rounding is consistent
for icond=1:length(events)
    events(icond).samples = round(events(icond).times ./ dt);
    events(icond).times   = events(icond).samples .* dt;
end

%% Build a design matrix to get basic signals
hrf_types = process_nst_glm_fit('get_hrf_types');
[X,names] = process_nst_glm_fit('make_design_matrix', time, events, ...
                                hrf_types.CANONICAL, 25, 0);
assert(all(strcmp(names, {events.label})));

%% Forge nirs data file
nirs.ml = [1 1 1 1;... %S1D1HbO
           1 1 1 2;... %S1D1HbR
           2 2 1 1;... %S2D2HbO
           2 2 1 2;... %S2D2HbR
           3 3 1 1;... %S3D3HbO
           3 3 1 2];   %S3D3HbR

nirs.SD.Lambda = [690, 830]; % Use wavelengths here to enable head model computation
                             % Even though it's Hb concentration data.
                             % The head model will only be used to get cortical
                             % vertices that are in the field of view
nirs.SD.SrcPos = [ 0.1073 -0.0399 0.0300;  %S1 
                  -0.0843 -0.0254 0.0643;  %S2
                   0.0269  0.0236 0.1402]; %S3

nirs.SD.DetPos = [ 0.1049 -0.0378 0.0585;  %D1 
                  -0.0863  0.0031 0.0681;  %D2
                   0.0513  0.0106 0.1348]; %D3
            
           % stim1   stim2  
delta_hb = [ 20   15     ;... %HbO (mumol.l-1)
            -15  -10    ];    %HbR (mumol.l-1)

nirs.d = randn(length(time), size(nirs.ml, 1)) * 0.5;
nirs.d(:, 1) = nirs.d(:, 1) + X(:, 1) .* delta_hb(1,1);
nirs.d(:, 2) = nirs.d(:, 2) + X(:, 1) .* delta_hb(2,1);

nirs.d(:, 3) = nirs.d(:, 3) + X(:, 2) .* delta_hb(1,2);
nirs.d(:, 4) = nirs.d(:, 4) + X(:, 2) .* delta_hb(2,2);
nirs.t = time;

save(nirs_data_fn, '-struct', 'nirs', '-mat');

% Create subject with Colin27 anantomy
[subject_name, sSubject, iSubject] = bst_create_test_subject('Colin27_4NIRS_lowres', 0);

% Import nirs data
sNirs = utest_import_nirs_in_bst(nirs_data_fn, 0, 0);

% Rename to Hb
sNirs = bst_process('CallProcess', 'process_set_comment', sNirs, [], ...
                    'tag', 'Hb', 'isindex', 0);

% Inject events
process_nst_import_csv_events('import_events', [], sNirs, events);

% Compute headmodel
bst_process('CallProcess', 'process_nst_import_head_model', sNirs, [], ...
            'use_closest_wl', 1, 'use_all_pairs', 0, ...
            'force_median_spread', 0, ...
            'normalize_fluence', 0, ...
            'sensitivity_threshold_pct', 0.5, ...
            'smoothing_fwhm', 0);

% Set unit to mumol
sData = load(file_fullpath(sNirs.FileName));
sData.DisplayUnits = 'mumol.l-1';
save(file_fullpath(sNirs.FileName), '-struct', 'sData');        

% Replace WL with Hb
ChannelMat = in_bst_channel(sNirs.ChannelFile);
ChannelMat.Nirs = rmfield(ChannelMat.Nirs, 'Wavelengths');
ChannelMat.Nirs.Hb = {'HbO', 'HbR'};
[tmp, iChannelStudy] = bst_get('ChannelForStudy', sNirs.iStudy);
db_set_channel(iChannelStudy, ChannelMat, 0, 0);
    
% Create cortical signals (duplicate channel signals)
sStudy = bst_get('Study', sNirs.iStudy);
head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
vertices_S1D1 = squeeze(head_model.Gain(strcmp(head_model.pair_names, 'S1D1'), 1, :) > 0);
vertices_S2D2 = squeeze(head_model.Gain(strcmp(head_model.pair_names, 'S2D2'), 1, :) > 0);
vertices_S3D3 = squeeze(head_model.Gain(strcmp(head_model.pair_names, 'S3D3'), 1, :) > 0);


nb_vertices = size(head_model.Gain, 3);
nb_samples = length(time);
pdata_hbo = zeros(nb_vertices, nb_samples);
pdata_hbo(vertices_S1D1, :) = repmat(nirs.d(:, 1)', sum(vertices_S1D1), 1);
pdata_hbo(vertices_S2D2, :) = repmat(nirs.d(:, 3)', sum(vertices_S2D2), 1);
pdata_hbo(vertices_S3D3, :) = repmat(nirs.d(:, 5)', sum(vertices_S3D3), 1);

pdata_hbr = zeros(nb_vertices, nb_samples);
pdata_hbr(vertices_S1D1, :) = repmat(nirs.d(:, 2)', sum(vertices_S1D1), 1);
pdata_hbr(vertices_S2D2, :) = repmat(nirs.d(:, 4)', sum(vertices_S2D2), 1);
pdata_hbr(vertices_S3D3, :) = repmat(nirs.d(:, 5)', sum(vertices_S3D3), 1);

extra.DisplayUnits = 'mumol.l-1';
nst_bst_add_surf_data(pdata_hbo, sData.Time, ...
    [], 'hbo_cortex', 'HbO cortex', ...
    [], sStudy,  'Cortical HbO signals', head_model.SurfaceFile, 1, extra);

nst_bst_add_surf_data(pdata_hbo, sData.Time, ...
    [], 'hbr_cortex', 'HbR cortex', ...
    [], sStudy,  'Cortical HbR signals', head_model.SurfaceFile, 1, extra);


% Remove head model
delete_head_model(sStudy, sNirs.iStudy, 1);

panel_protocols('ReloadNode', panel_protocols('SelectStudyNode', sNirs.iStudy));
db_save();

% Save all subject data
subject_data_fn = fullfile(dest_dir, 'activity_test_subject.zip');
export_protocol(bst_get('iProtocol'), iSubject, subject_data_fn);

% Save ground-truth activity info
activity_info(1).condition_name = stim_event_names{1};
activity_info(1).activ_chans = [1 2];
activity_info(1).activ_scout.Handles = [];
activity_info(1).activ_scout.Region = 'RF';
activity_info(1).activ_scout.Function = 'mean';
activity_info(1).activ_scout.Label = 'activ_stim1';
activity_info(1).activ_scout.Color = [1 0 0];
activity_info(1).activ_scout.Seed = 4948;
activity_info(1).activ_scout.Vertices = [4905 4912 4921 4922 4941 4942 4948 4949 4959 4960 4961 4965 4973 4977];

activity_info(2).condition_name = stim_event_names{2};
activity_info(2).activ_chans = [3 4];
activity_info(2).activ_scout.Handles = [];
activity_info(2).activ_scout.Region = 'RO';
activity_info(2).activ_scout.Function = 'mean';
activity_info(2).activ_scout.Label = 'activ_stim2';
activity_info(2).activ_scout.Color = [0 0.8 0];
activity_info(2).activ_scout.Seed = 2514;
activity_info(2).activ_scout.Vertices = [2506 2508 2509 2511 2514 2515 2516 2520 2521 2525 2526 2530 2533 2539 2545 2552 2553 2562 2571];

activity_info_fn = fullfile(dest_dir, 'activity_info.mat');
save(activity_info_fn, 'activity_info', '-mat');

end


function sStudy = delete_head_model(sStudy, iStudy, iHeadModelDel)

file_delete(file_fullpath(sStudy.HeadModel(iHeadModelDel).FileName), 1);

% From node_delete.m / case 'headmodel'

% Remove files descriptions from database
sStudy.HeadModel(iHeadModelDel) = [];
% Update default headmodel
nbHeadModel = length(sStudy.HeadModel);
if (nbHeadModel <= 0)
    sStudy.iHeadModel = [];
elseif (nbHeadModel == 1)
    sStudy.iHeadModel = 1;
elseif (sStudy.iHeadModel > nbHeadModel)
    sStudy.iHeadModel = nbHeadModel;
else
    % Do not change iHeadModel
end
% Study was modified
bst_set('Study', iStudy, sStudy);
panel_protocols('UpdateNode', 'Study', iStudy);
db_save();
end