function OutputFile = bst_create_nirs_data(condition_label, signals, time, chan_names, srcs_pos, dets_pos)
% BST_CREATE_NIRS_DATA forge NIRS data into brainstorm for the subject "test_subject" in the protocol "nst_utest".
%
%   [OutputFile] = BST_CREATE_NIRS_DATA(CONDITION_LABEL, SIGNALS, TIME, CHAN_NAMES, SRCS_POS, DETS_POS)
%       When a parameter is not given, the following defaults are used:
%           - 4 channels: 'S1D1WL685', 'S1D1WL830', 'S3D4WL685', 'S3D4WL830'
%           - 1000 time points with a temporal resolution of 10 Hz (-> 100 sec)
%       CONDITION_LABEL (str): name of the condition where to put the forged data.
%          If not given or empty, default is 'dummy_nirs_data'
%       SIGNALS (matrix of double): nirs signals, size: (nb_channels x nb_samples)
%          If not given or empty, default is randn(nb_channels x nb_samples) 
%       TIME (array of double): temporal positions in sec.
%          If not given or empty, default is 0:0.1:nb_samples
%       CHAN_NAMES (cellarray of str): each channel formatted as 'SxDyWLz' or 'SxDyHbt'
%           where:
%               x: source index
%               y: detector index
%               z: wavelength
%               t: Hb type (O, R, T).
%           Examples: S1D2WL685, S01D7WL830, S3D01HbR
%           If not given of empty, default is {'S1D1WL685', 'S1D1WL830','S3D4WL685', 'S3D4WL830'}
%       SRCS_POS (matrix of double): coordinates of channel sources, size: (3 x nb_channels)
%           Note: aligns with ChannelMat.Channel(:).Loc(:,1)
%       DETS_POS (matrix of double): coordinates of channel detectors, size: (3 x nb_channels)
%           Note: aligns with ChannelMat.Channel(:).Loc(:,2)
%
%       OUTPUTFILE(struct): brainstorm data outputFile, see db_template('data').
%
%   See also NST_FORMAT_CHANNELS, NST_FORMAT_CHANNEL

DEFAULT_WAVELENGTHS = [685 830]; %nm
DEFAULT_DT = 0.1; %sec
DEFAULT_NB_SAMPLES = 1000;
DEFAULT_CONDITION = 'dummy_nirs_data';
DEFAULT_CHANNELS = {'S1D1WL685', 'S1D1WL830', 'S3D4WL685', 'S3D4WL830'};

ProtocolName = 'nst_utest';
sProtocol = bst_get('Protocol', ProtocolName);
if isempty(sProtocol)
    error('Protocol "nst_utest" not found.');
end
gui_brainstorm('SetCurrentProtocol', sProtocol);

if nargin < 1 || isempty(condition_label)
   condition_label = DEFAULT_CONDITION; 
end

if nargin < 2 || isempty(signals)
   if nargin > 4 && ~isempty(srcs_pos)
       if nargin < 6 || isempty(dets_pos)
           error('"dets_pos" must be defined when "srcs_pos" is given');
       end
       nb_channels = size(src_pos, 2);
   else
       nb_channels = length(DEFAULT_CHANNELS);
       chan_names = DEFAULT_CHANNELS;
   end
   
   if nargin > 2 && ~isempty(time)
       nb_samples = length(time);
   else
       nb_samples = DEFAULT_NB_SAMPLES;
   end
   
   signals = randn(nb_channels, nb_samples);
end

if nargin < 3 || isempty(time)
    time = 0:DEFAULT_DT:size(signals, 2);
end

iStudy = db_add_condition('test_subject', condition_label);
sStudy = bst_get('Study', iStudy);

nb_channels = size(signals,1);

if nargin < 4 || isempty(chan_names)
    chan_names = cell(1, nb_channels);
    wavelengths = DEFAULT_WAVELENGTHS;
    for ichan=1:nb_channels
        i_src = floor((ichan-1) / 4) + 1;
        i_det = floor((ichan-1) / 2) + 1;
        i_wl = mod(ichan-1,2) + 1;
        chan_names{ichan} = nst_format_channel(i_src, i_det, wavelengths(i_wl));
    end
end
assert(length(chan_names)==nb_channels);

if nargin < 5 || isempty(srcs_pos)
    if nargin > 5 && ~isempty(dets_pos)
        error('"srcs_pos" must be defined when "dets_pos" is given');
        return;
    end
    [i_srcs, i_dets, measures, ctype] = nst_unformat_channels(chan_names);
    % Default locations:
    %   - one source every 3cm on a straight line at y=0, z=0
    %   - one detector every 3cm on a straight line at y=3, z=0
    %   - use an origin point at the back of the head 
    %    -> using (0,0,0) crashes brainstorm. This point does make sense
    %       anyway since it's in the middle of the head in the SCS referential
    srcs_pos = zeros(3, nb_channels);
    dets_pos = zeros(3, nb_channels);
    ref_point =  [-0.0789 ; -0.0263 ; 0.0655]; % meter, somewhere at the back of the head when using Colin27_4NIRS
    for ichan=1:nb_channels
       srcs_pos(:, ichan) = ref_point + [(i_srcs(ichan)-1) * 0.03; 0; 0]; %meter
       dets_pos(:, ichan) = ref_point + [(i_dets(ichan)-1) * 0.03; 0.03; 0]; %meter
    end
end

if nargin==5
    error('"dets_pos" must be defined when "srcs_pos" is given');
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