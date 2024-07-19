function  tutorial_nirstorn_2024(tutorial_dir, reports_dir)
% TUTORIAL_NIRSTORM_2024: Script that reproduces the results of the
% nirstorm article
%
% CORRESPONDING ONLINE TUTORIALS:
%     https://neuroimage.usc.edu/brainstorm/Tutorials/NIRSTORM
%
% INPUTS: 
%    - tutorial_dir: Directory where the sample_nirstorm.zip file has been downloaded
%    - reports_dir  : Directory where to save the execution report (instead of displaying it)

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Author: Edouard Delaire, 2024

% ===== FILES TO IMPORT =====
% Output folder for reports
if (nargin < 2) || isempty(reports_dir) || ~isdir(reports_dir)
    reports_dir = [];
end
% You have to specify the folder in which the tutorial dataset is unzipped
if (nargin == 0) || isempty(tutorial_dir) || ~file_exist(tutorial_dir)
    error('The first argument must be the full path to the tutorial dataset folder.');
end

% Processing options: 
filter_type     = 'FIR'; % Options : {'FIR', 'IIR', 'NONE'}
filter_band     = [0.005, 0.1];
epoch_duration  = [-10, 35] ; 


SubjectName = 'sub-01';

participants_file           = fullfile(tutorial_dir, 'participants.tsv');
FidFile                     = fullfile(tutorial_dir, SubjectName, 'anat','T1w.json');
FSPath                      = fullfile(tutorial_dir, 'derivatives', 'FreeSurfer/',SubjectName);
TissueSegmentationFile      = fullfile(tutorial_dir, 'derivatives', 'segmentation',SubjectName,'segmentation_5tissues.nii');
scoutPath                   = fullfile(tutorial_dir, 'derivatives', 'segmentation',SubjectName,'scout_hand.label');

RawFile     = fullfile(tutorial_dir, SubjectName, 'nirs',sprintf('%s_task-tapping_run-01.snirf',SubjectName));
FluenceDir  = fullfile(tutorial_dir, 'derivatives','Fluences');

% Check if the folder contains the required files
% if ~file_exist(RawFile)
%     error(['The folder ' tutorial_dir ' does not contain the folder from the file sample_epilepsy.zip.']);
% end

sParticipant = readtable(participants_file, 'FileType','text');

% ===== CREATE PROTOCOL =====
ProtocolName = 'TutorialNIRSTORM';

% Start brainstorm without the GUI
if ~brainstorm('status')
    brainstorm nogui
end

% Delete existing protocol
%gui_brainstorm('DeleteProtocol', ProtocolName);
% Create new protocol
%gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

% Start a new report
bst_report('Start');

%% Part 1. Importing subject anatomical data
json = bst_jsondecode(FidFile);
sFid = process_import_bids('GetFiducials',json, 'voxel');

    bst_process('CallProcess', 'process_import_anatomy', [], [], ...
        'subjectname', SubjectName, ...
        'mrifile',     {FSPath, 'FreeSurfer+Thick'}, ...
        'nvertices',   15000, ...
        'nas',         sFid.NAS, ...
        'lpa',         sFid.LPA, ...
        'rpa',         sFid.RPA, ...
        'ac',          sFid.AC, ...
        'pc',          sFid.PC, ...
        'ih',          sFid.IH);

[sSubject, iSubject]  = bst_get('Subject', SubjectName);

% Import tissue segmentation
tissues =  mri_getlabels('tissues5');
tissues(6) = {1}; tissues(5) = {2}; tissues(4) = {3}; tissues(3) = {4}; tissues(2) = {5};  
import_mri(iSubject, TissueSegmentationFile, '', 0, 0, 'segmentation_5tissues',tissues);

% Remesh skin to 10 000 vertices, iso-mesh
tess_remesh(sSubject.Surface(sSubject.iScalp).FileName ,10000 )

% Choose mid-surface as default cortical surface
db_surface_default(iSubject, 'Cortex', find(contains({sSubject.Surface.FileName},'tess_cortex_mid_low.mat')));
panel_protocols('RepaintTree');

% Import hand-know region
import_label(sSubject.Surface(sSubject.iCortex).FileName, scoutPath,0);

% Compute voronoi-interpolation
bst_process('CallProcess', 'process_nst_compute_voronoi', [], [], ...
    'subjectname',  SubjectName, ...
    'do_grey_mask', 1);



%% Part 2. Import functional data


sFile = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
                                    'subjectname',    SubjectName, ...
                                    'datafile',       {RawFile, 'NIRS-SNIRF'}, ...
                                    'channelreplace', 1, ...
                                    'channelalign',   1, ...
                                    'evtmode',        'value');


% Import head points 
% Note: already done since headpoints are included in the snirf file
% bst_process('CallProcess', 'process_headpoints_add', sFile, [], ...
%             'channelfile', {HeadpointsPath, 'ASCII_NXYZ'}, ...
%             'fixunits',    0.1, ...
%             'vox2ras',     0); % we are using fiducials define in headpoints 

% Remove headpoints bellow nasion
% Note: already done since headpoints are included in the snirf file
% process_headpoints_remove('RemoveHeadpoints', sFile.ChannelFile, 0) 

% Refine registration
% Note: already done since headpoints are included in the snirf file
% bst_process('CallProcess', 'process_headpoints_refine', sFile, []);

% Process: Snapshot: Sensors/MRI registration
bst_process('CallProcess', 'process_snapshot', sFile, [], ...
    'target',   1, ...  % Sensors/MRI registration
    'modality', 7, ...  % NIRS
    'orient',   1, ...  % left
    'Comment',  'NIRS/MRI Registration');

% Process: Duplicate tapping events
sFile = bst_process('CallProcess', 'process_evt_merge', sFile, [], ...
    'evtnames', 'tapping', ...
    'newname',  'tapping/start', ...
    'delete',   0);

% Process: Convert to simple event
sFile = bst_process('CallProcess', 'process_evt_simple', sFile, [], ...
    'eventname',  'tapping/start', ...
    'method', 'start');

% Process: Detect bad channels
sFile = bst_process('CallProcess', 'process_nst_detect_bad', sFile, [], ...
    'option_sci',                   0, ...
    'sci_threshold',                80, ...
    'power_threshold',              10, ...
    'option_coefficient_variation', 1, ...
    'coefficient_variation',        10, ...
    'option_remove_saturating',     0, ...
    'option_max_sat_prop',          10, ...
    'option_min_sat_prop',          10, ...
    'option_separation_filtering',  0, ...
    'option_separation',            [0, 5], ...
    'auxilary_signal',              3, ...  % Remove all
    'option_keep_unpaired',         0);

sRaw = sFile;

%% Part 3. Preprocessing

sRawdOD = bst_process('CallProcess', 'process_nst_dOD', sRaw, [], ...
    'option_baseline_method', 1, ...  % mean
    'timewindow',             []);

sRawdODFiltered = bst_process('CallProcess', 'process_bandpass', sRawdOD, [], ...
        'sensortypes', 'NIRS', ...
        'highpass',    filter_band(1), ...
        'lowpass',     filter_band(2), ...
        'tranband',    0.005, ...
        'attenuation', 'relax', ...     % 40dB (relaxed)
        'ver',         '2019', ...      % 2019
        'mirror',      0, ...
        'overwrite',   0);

sRawdODFilteredSS = bst_process('CallProcess', 'process_nst_remove_ssc', sRawdODFiltered, [], ...
    'SS_chan',                 'name', ...  % Based on Names
    'SS_chan_name',            'S1D17,S2D17', ...
    'separation_threshold_cm', 1.5);


%% Part 4.  Estimation of the HRF to tapping 

sTrialsOd = bst_process('CallProcess', 'process_import_data_event', sRawdODFilteredSS, [], ...
    'subjectname', SubjectName, ...
    'condition',   '', ...
    'eventname',   'tapping/start', ...
    'timewindow',  [], ...
    'epochtime',   epoch_duration, ...
    'createcond',  0, ...
    'ignoreshort', 0, ...
    'usectfcomp',  0, ...
    'usessp',      0, ...
    'freq',        [], ...
    'baseline',    []);

trialStatus = true(1,20);
trialStatus([1, 15]) = false;

% Process: DC offset correction: [-10.000s,0.000s]
sTrialsOd = bst_process('CallProcess', 'process_baseline_norm', sTrialsOd(trialStatus), [], ...
    'baseline',    [-10, 0], ...
    'sensortypes', 'NIRS', ...
    'method',      'bl', ...  % DC offset correction:    x_std = x - &mu;
    'overwrite',   0);

% Process: Average+Stderr: By trial group (folder average)
sAverageOd = bst_process('CallProcess', 'process_average', sTrialsOd, [], ...
    'avgtype',       5, ...  % By trial group (folder average)
    'avg_func',      7, ...  % Arithmetic average + Standard error
    'weighted',      0, ...
    'keepevents',    1);

sAverageHb  = bst_process('CallProcess', 'process_nst_mbll_dOD', sAverageOd, [], ...
    'option_age',         sParticipant{1,'age'}, ...
    'option_pvf',         50, ...
    'option_do_plp_corr', 1, ...
    'option_dpf_method',  1);  % SCHOLKMANN2013

%% Part 5.  Localization of the response on the cortex 

bst_process('CallProcess', 'process_nst_cpt_fluences', sFilesOD, [], ...
    'subjectname',  SubjectName, ...
    'fluencesCond', struct(...
         'surface',                   'montage', ...
         'ChannelFile',               sFilesOD.ChannelFile    , ...
         'SubjectName',               SubjectName, ...
         'segmentation_label',        2, ...
         'wavelengths',               '685 ,830', ...
         'software',                  'mcxlab-cuda', ...
         'mcxlab_gpuid',              1, ...
         'mcxlab_nphoton',            100, ...
         'outputdir',                 FluenceDir, ...
         'mcxlab_flag_thresh',        0, ...
         'mcxlab_overwrite_fluences', 0, ...
         'mcxlab_flag_autoOP',        1));

% Save and display report
ReportFile = bst_report('Save', []);
if ~isempty(reports_dir) && ~isempty(ReportFile)
    bst_report('Export', ReportFile, reports_dir);
else
    bst_report('Open', ReportFile);
end

disp([10 'BST> tutorial_nirstorm: Done.' 10]);

end

