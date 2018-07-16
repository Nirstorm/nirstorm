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
            
            dt = 0.1; %sec (10Hz)
            time = (0:100) * dt; %sec
            nb_samples = length(time);
            
            bst_create_test_subject()
            
            %% Only noise -- uncorrelated wavelength signals
            signals = [randn(1, nb_samples);randn(1, nb_samples)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            outputs = bst_process('CallProcess', 'process_nst_sci', [], nirs_input);
            %Check sci==0
            
            %% Perfect cardiac consistency
            
            %Check sci==1
            
            
            %% Noisy consistent cardiac components - high SNR
            
            %Check sci > 0.75
            
            %% Noisy consistent cardiac components - very low SNR
            
            %Check sci < 0.75
            
            
        end
        
        function test_sci_on_tapping_data(testCase)
            
        end
        
        
    end
end

function OutputFile = bst_create_nirs_data(label, signals, time, chan_names, srcs_pos, dets_pos)
% TODO: expose, doc, test
iStudy = db_add_condition('nst_utest', 'test_sci');
sStudy = bst_get('Study', iStudy);

nb_channels = size(signals,1);

if nargin < 4 || isempty(chan_names)
    chan_names = cell(1, nb_channels);
    for ichan=1:nb_channels
        i_src = round((ichan-1) / 4) + 1;
        i_det = round((ichan-1) / 2) + 1;
        i_wl = mod(ichan-1,2) + 1;
        chan_names{ichan} = nst_format_channel(i_src, i_det, i_wl);
    end
end
assert(length(chan_names)==nb_channels);

if nargin < 5 || isempty(srcs_pos)
    [i_srcs, i_dets, i_wls] = nst_unformat_channels(chan_names); %TODO - also make some checks for consistency
    % Default locations: 
    %   one source every 3cm on a straight line at y=0, z=0
    %   one detector every 3cm on a straight line at y=3, z=0
    for ichan=1:nb_channels
       srcs_pos = [i_srcs(ichan)-1 * 0.03; 0; 0]; %meter
       dets_pos = [i_dets(ichan)-1 * 0.03; 0.03; 0]; %meter
    end
end

if nargin < 6
    error('Provide det_pos with src_pos');
end

% Create channel definition
ChannelMat = db_template('channelmat');
ChannelMat.Channel = repmat(db_template('channeldesc'), 1, nb_channels);
for ichan=1:nb_channels
    [isrc, idet] = unformat_channel_name(chan_names{ichan});
    ChannelMat.Channel(ichan).Name = chan_names{ichan};
    ChannelMat.Channel(ichan).Group = 'NIRS';
    ChannelMat.Channel(ichan).Loc(:,1) = srcs_pos(:,ichan);
    ChannelMat.Channel(ichan).Loc(:,2) = dets_pos(:,ichan);
end

% Save channel definition
[tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
db_set_channel(iChannelStudy, ChannelMat, 0, 0);

% Save time-series data
sDataOut = db_template('data');
sDataOut.F            = signals;
sDataOut.Comment      = label;
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
end