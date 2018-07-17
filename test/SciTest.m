classdef SciTest < matlab.unittest.TestCase
    
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
        
        function test_sci_with_simulations(testCase)
            
            dt = 0.05; %sec (20Hz)
            time = (0:8000) * dt; %sec
            nb_samples = length(time);
            
            bst_create_test_subject();
            
            %% Only noise -- uncorrelated wavelength signals
            signals = [randn(1, nb_samples);randn(1, nb_samples);];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI close to 0
            testCase.assertTrue(all(sDataOut.F<0.1));
            
            %% Correlation outside of cardiac band
            freq = 5; % Hz
            signals = [cos(4*pi*freq*time) + randn(1, nb_samples);...
                       cos(4*pi*freq*time) + randn(1, nb_samples)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI close to 0
            testCase.assertTrue(all(sDataOut.F<0.1));
            
            %% Perfect cardiac consistency - correlated
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time); cos(4*pi*freq*time)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI high
            testCase.assertTrue(all(sDataOut.F>0.99));
            
            %% Perfect cardiac consistency - anticorrelated
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time); -cos(4*pi*freq*time)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI high
            testCase.assertTrue(all(sDataOut.F>0.99));
                        
            %% Consistent cardiac components - high SNR
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time) + randn(1, nb_samples)*0.2; ...
                       -cos(4*pi*freq*time) + randn(1, nb_samples)*0.2];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI > 0.75
            testCase.assertTrue(all(sDataOut.F>0.75));
            
            %% Consistent cardiac components - very low SNR
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time) + randn(1, nb_samples)*50; ...
                       -cos(4*pi*freq*time) + randn(1, nb_samples)*50];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI < 0.75
            testCase.assertTrue(all(sDataOut.F<0.75));
            
            
        end
        
        function test_sci_on_tapping_data(testCase)
            repo_url = nst_get_repository_url();
            data_fns = nst_request_files({{'unittest','motor_data','motor.nirs'}, ...
                                          {'unittest','motor_data','optodes.txt'}, ...
                                          {'unittest','motor_data','fiducials.txt'},...
                                          {'unittest','motor_data','sci.mat'}}, ...
                                         1, repo_url);

                                    
            nirs_fn = data_fns{1};
            sFile = utest_import_nirs_in_bst(nirs_fn);
            output = bst_process('CallProcess', 'process_nst_sci', sFile, []);
            sDataOut = in_bst_data(output.FileName);
            
            % Non-regression test (sci results manually checked and
            % recorded on 17 July 2018)
            expected_sci_fn = data_fns{end};
            expected_sci = load(expected_sci_fn);
            expected_sci = expected_sci.sci_motor;
            
            testCase.assertTrue(all_close(sDataOut.F, expected_sci, 0.01, 0.01));
        end
        
        
    end
end

function OutputFile = bst_create_nirs_data(label, signals, time, chan_names, srcs_pos, dets_pos)
% TODO: expose, doc, test
iStudy = db_add_condition('test_subject', label);
sStudy = bst_get('Study', iStudy);

nb_channels = size(signals,1);

if nargin < 4 || isempty(chan_names)
    chan_names = cell(1, nb_channels);
    wavelengths = [685 830];
    for ichan=1:nb_channels
        i_src = floor((ichan-1) / 4) + 1;
        i_det = floor((ichan-1) / 2) + 1;
        i_wl = mod(ichan-1,2) + 1;
        chan_names{ichan} = nst_format_channel(i_src, i_det, wavelengths(i_wl));
    end
end
assert(length(chan_names)==nb_channels);

if nargin < 5 || isempty(srcs_pos)
    [i_srcs, i_dets, measures, ctype] = nst_unformat_channels(chan_names);
    % Default locations: 
    %   one source every 3cm on a straight line at y=0, z=0
    %   one detector every 3cm on a straight line at y=3, z=0
    srcs_pos = zeros(3, nb_channels);
    dets_pos = zeros(3, nb_channels);
    ref_point =  [-0.0789 ; -0.0263 ; 0.0655]; % meter, somewhere at the back of the head when using Colin27_4NIRS
    for ichan=1:nb_channels
       srcs_pos(:, ichan) = ref_point + [(i_srcs(ichan)-1) * 0.03; 0; 0]; %meter
       dets_pos(:, ichan) = ref_point + [(i_dets(ichan)-1) * 0.03; 0.03; 0]; %meter
    end
end

if nargin>=5 && nargin < 6
    error('Provide det_pos with src_pos');
end

% Create channel definition
[i_srcs, i_dets, measures, ctype] = nst_unformat_channels(chan_names);
chan_types = nst_channel_types();
ChannelMat = db_template('channelmat');
ChannelMat.Comment = sprintf('NIRS-BRS channels (%d)', nb_channels);
ChannelMat.Channel = repmat(db_template('channeldesc'), 1, nb_channels);
for ichan=1:nb_channels
    ChannelMat.Channel(ichan).Name = chan_names{ichan};
    ChannelMat.Channel(ichan).Type = 'NIRS';
    ChannelMat.Channel(ichan).Weight = 1;
    if ctype == chan_types.WAVELENGTH
        ChannelMat.Channel(ichan).Group = sprintf('WL%d', measures(ichan));
    elseif ctype == chan_types.HB
        ChannelMat.Channel(ichan).Group = measures{ichan};
    end
    ChannelMat.Channel(ichan).Loc(:,1) = srcs_pos(:, ichan);
    ChannelMat.Channel(ichan).Loc(:,2) = dets_pos(:, ichan);
end
if ctype == chan_types.WAVELENGTH
    ChannelMat.Nirs.Wavelengths = unique(measures);
elseif ctype == chan_types.HB
    ChannelMat.Nirs.Hb = unique(measures);
end

% Save channel definition
[tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
db_set_channel(iChannelStudy, ChannelMat, 0, 0);

% Save time-series data
sDataOut = db_template('data');
sDataOut.F            = signals;
sDataOut.Comment      = 'data';
sDataOut.ChannelFlag  = ones(nb_channels, 1);
sDataOut.Time         = time;
sDataOut.DataType     = 'recordings';
sDataOut.nAvg         = 1;
sDataOut.DisplayUnits = 'mol.l-1'; %TODO: expose this?

% Generate a new file name in the same folder
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_forged');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(iStudy, OutputFile, sDataOut);
db_save();
end


function flag = all_close(v1, v2, rtol, atol)

if nargin < 3
    rtol = 1e-5; %default relative tolerance
end

if nargin < 4
    atol = 1e-5; %default absolute tolerance
end

% Convert to vectors if needed:
if ndims(v1) == 2
    v1 = v1(:)';
end

if ndims(v2) == 2
    v2 = v2(:)';
end

flag = all(abs(v1 - v2) <= (atol + rtol * max(abs(v1), abs(v2))));

end