function nst_tutorial_tapping(tutorial_dir, reports_dir)
% TUTORIAL_INTRODUCTION: Script that runs the tutorial on basic NIRS data
% processing (finger tapping experiment). 
% See http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping
%
% INPUTS: 
%    - tutorial_dir: Directory where the sample_nirs.zip file has been unzipped
%
% Author: Thomas Vincent, 2017

% Output folder for reports
if (nargin < 2) || isempty(reports_dir) || ~isdir(reports_dir)
    reports_dir = [];
end
% You have to specify the folder in which the tutorial dataset is unzipped
if (nargin == 0) || isempty(tutorial_dir) || ~file_exist(tutorial_dir)
    error('The first argument must be the full path to the dataset folder.');
end

SubjectName = 'Subject01';
% Build the path of the files to import
AnatDir    = fullfile(tutorial_dir, 'sample_nirs', 'anatomy');
Run1File   = fullfile(tutorial_dir, 'sample_nirs', 'data', ...
                      'S01_Block_FO_LH_Run01.nirs');
% Check if the folder contains the required files
if ~file_exist(Run1File)
    error(['The folder ' tutorial_dir ' does not contain the folder from the file sample_nirs.zip.']);
end

disp([10 'DEMO> Tutorial NIRSFingerTapping: Create protocol' 10]);
ProtocolName = 'TutorialNIRS';
% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end

% Delete existing protocol
gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);
% Start a new report
bst_report('Start');

disp([10 'DEMO> Tutorial NIRSFingerTapping: Import anatomy' 10]);
bst_process('CallProcess', 'process_import_anatomy', [], [], ...
    'subjectname', SubjectName, ...
    'mrifile',     {AnatDir, 'BrainVISA'}, ...
    'nvertices',   15000, ...
    'nas', [95 213 114], ...
    'lpa', [31 126 88 ], ...
    'rpa', [164 128 89]);

disp([10 'DEMO> Tutorial NIRSFingerTapping: Import data' 10]);
% Process: Create link to raw files
sFilesRun1 = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
    'subjectname',  SubjectName, ...
    'datafile',     {Run1File, 'NIRS-BRS'}, ...
    'channelalign', 1);
% Process: Read from channel
bst_process('CallProcess', 'process_evt_read', sFilesRun1, [], ...
    'stimchan',  'NIRS_AUX', ...
    'trackmode', 3, ...  % Value: detect the changes of channel value
    'zero',      0);
% Process: Rename event (standard_fix>standard)
bst_process('CallProcess', 'process_evt_rename', sFilesRun1, [], ...
    'src',  'AUX1', ...
    'dest', 'MOTOR');

% Define movement events
movement_events = db_template('event');
movement_label = 'NIRS_mvt';
movement_events.label = movement_label;
movement_events.times = [577.3 662.2 689.3; ...
                         580.4 663.6 691.6];
movement_events.epochs = ones(1, size(movement_events.times, 2));
%TODO: Use more generic process?
process_nst_import_csv_events('import_events', [], sFilesRun1, movement_events);

sFileMoco = bst_process('CallProcess', 'process_nst_motion_correction', sFilesRun1, [], ...
    'option_event_name', movement_label);

bst_process('CallProcess', 'process_nst_detect_bad', sFileMoco, [], ...
            'option_remove_negative', 1, ...
            'option_invalidate_paired_channels', 1, ...
            'option_max_sat_prop', 1);

sFileHb = bst_process('CallProcess', 'process_nst_mbll', sFileMoco, [], ...
                        'option_age', {25, ''}, ...
                        'option_pvf', {50, ''}, ...
                        'option_baseline_method', {1, {'mean', 'median'}}, ...
                        'option_do_plp_corr', 1);

sFileHb_detrend = bst_process('CallProcess', 'process_detrend', sFileHb, [], ...
                              'timewindow', [], ...
                              'sensortypes', 'NIRS');
        
sFileHb_bp = bst_process('CallProcess', 'process_bandpass', sFileHb_detrend, [], ...
                              'highpass', {0.005, ''}, ...
                              'lowpass', {0.08, ''}, ...
                              'attenuation', 'relax', ...
                              'mirror', 0, ...
                              'sensortypes', 'NIRS');

%% Import peristimuli blocks
% Read channel file
ChannelFile = bst_get('ChannelFileForStudy', sFileHb_bp.FileName);
ChannelMat = in_bst_channel(ChannelFile);
% Get sFile structure
sFile = in_fopen(sFileHb_bp.FileName, 'BST-DATA');

ImportOptions = db_template('ImportOptions');
ImportOptions.ImportMode = 'Event';
ImportOptions.UseEvents = 1;
ImportOptions.TimeRange = sFile.prop.times;
ImportOptions.EventsTimeRange = [-5.0, 55.0]; %sec
ImportOptions.GetAllEpochs = 0;
ImportOptions.iEpochs = 1;
ImportOptions.SplitRaw = 0;
ImportOptions.SplitLength = [];
ImportOptions.Resample = 0;
ImportOptions.ResampleFreq = [];
ImportOptions.events = sFile.events(strcmp({sFile.events.label}, 'MOTOR')); 
ImportOptions.UseCtfComp = 1;          
ImportOptions.UseSsp = 0;
ImportOptions.RemoveBaseline = 'time';
ImportOptions.BaselineRange = [-5.0, -0.1]; %sec
ImportOptions.IgnoreShortEpochs = 0;
ImportOptions.DisplayMessages = 0;
[sSubject, iSubject] = bst_get('Subject', SubjectName);
sMotorBlocks = import_data(sFile, ChannelMat, sFile.format, [], iSubject, ImportOptions);

% Compute average with std dev.
sFilesAvg = bst_process('CallProcess', 'process_average', sMotorBlocks, [], ...
    'avgtype',    1, ...  % All files
    'avg_func',   6, ...  % avg + std dev
    'weighted',   0, ...
    'keepevents', 0);
hFigNirs1  = view_timeseries(sFilesAvg(1).FileName, 'NIRS');
hFigTopo1 = view_topography(sFilesAvg(1).FileName, 'NIRS', '3DOptodes');
panel_time('SetCurrentTime', 13.0);

end