function  tutorial_nirstorn_2024(tutorial_dir, reports_dir)
% TUTORIAL_EPILEPSY: Script that reproduces the results of the online tutorial "EEG/Epilepsy".
%
% CORRESPONDING ONLINE TUTORIALS:
%     https://neuroimage.usc.edu/brainstorm/Tutorials/Epilepsy
%
% INPUTS: 
%    - tutorial_dir: Directory where the sample_epilepsy.zip file has been unzipped
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
% Author: Edouard Delaire, 2014-2022


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
filter_type = 'FIR'; % Options : {'FIR', 'IIR', 'NONE'}
filter_band = [0.005, 0.1];

epoch_duration = [-10, 35] ; 


ProtocolName = 'TutorialNIRSTORM';
SubjectName = 'sub-01';

FSPath      = fullfile(tutorial_dir, 'derivatives', 'FreeSurfer/','sub-01/');
RawFile    = fullfile(tutorial_dir, 'sub-01/', 'nirs','sub-01_task-tapping_run-01.snirf');
FidFile    = fullfile(tutorial_dir, 'sub-01/', 'anat','T1w.json');

scoutPath   = fullfile(tutorial_dir, 'derivatives', 'segmentation','sub-01','scout_hand.label');
TissueSegmentationFile   = fullfile(tutorial_dir, 'derivatives', 'segmentation','sub-01','segmentation_5tissues.nii');

MotionFile  = fullfile(tutorial_dir, 'workshop', 'Functional','events_tapping_motion.txt');
FluenceDir  = fullfile(tutorial_dir, 'derivatives');

% Check if the folder contains the required files
% if ~file_exist(RawFile)
%     error(['The folder ' tutorial_dir ' does not contain the folder from the file sample_epilepsy.zip.']);
% end


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
[sAllAtlas, err] = import_label(sSubject.Surface(sSubject.iCortex).FileName, scoutPath,0);

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


end

