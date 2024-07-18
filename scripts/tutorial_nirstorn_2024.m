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

ProtocolName = 'TutorialNIRSTORM';
SubjectName  = 'sub-01';

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

SubjectName = 'sub-01';

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





end

