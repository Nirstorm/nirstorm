function varargout = process_nst_compute_voronoi( varargin )

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Thomas Vincent, ZhengChen Cai (2017)
%          Edouard Delaire (2025)

eval(macro_method);
end

function sProcess = GetDescription() 
    % Description the process
    sProcess.Comment     = 'Compute Voronoi volume-to-cortex interpolator';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1401;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'import'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 0;
    
    % Option: Subject name
    sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = '';

    sProcess.options.do_grey_mask.Comment = 'Grey matter masking';
    sProcess.options.do_grey_mask.Type    = 'checkbox';
    sProcess.options.do_grey_mask.Value   = 1;      
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
OutputFiles = {};

% Get subject name
if isfield(sProcess.options, 'subjectname') && ~isempty(sProcess.options.subjectname.Value)
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
else
    SubjectName = sInputs.SubjectName;
end
if isempty(SubjectName)
    bst_report('Error', sProcess, [], 'Subject name is empty.');
    return;
end

% ===== GET SUBJECT =====
% Get subject
[sSubject, iSubject] = bst_get('Subject', SubjectName);
if isempty(iSubject)
    bst_report('Error', sProcess, [], ['Subject "' SubjectName '" does not exist.']);
    return
end

% Get given data information
voronoi_fn      = get_voronoi_fn(sSubject);
if exist(voronoi_fn,'file')
    if ~java_dialog('confirm', sprintf('File already exist.\nDo you want to overwrite it?'), 'Computing Voronoi...')
        return;
    end
end

seg_label       = 'segmentation_5tissues';
segmentation_id = find(strcmp(seg_label, {sSubject.Anatomy.Comment}));


bst_progress('start', 'MRI/Surface Voronoi interpolator','Computing Voronoi partitioning ...', 1, 2);

if ~isempty(segmentation_id) && sProcess.options.do_grey_mask.Value  


    sMri            = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
    sSegmentation   = in_mri_bst(sSubject.Anatomy(segmentation_id).FileName);

    if ~all(round(sMri.Voxsize(1:3) .* 1000) == round(sSegmentation.Voxsize(1:3) .* 1000))
        bst_report('Error', sProcess, [], 'MRI and Segmentation have different voxel size');
        return
    end

    [voronoi, sMRI] = Compute(sSubject.Surface(sSubject.iCortex).FileName, ...
                      sSubject.Anatomy(sSubject.iAnatomy).FileName, ...
                      sSubject.Anatomy(segmentation_id).FileName);

elseif isempty(segmentation_id) || ~sProcess.options.do_grey_mask.Value 
    
    msg = '';
    if sProcess.options.do_grey_mask.Value
        msg = sprintf('MRI segmentation (%s) not found. \n', seg_label);
    end
    msg = [ msg, 'Interpolator is not constrained to grey matter. Expect partial volume effect.'];
    
    bst_report('Warning', sProcess, sInputs, msg);

    [voronoi, sMRI] = Compute(sSubject.Surface(sSubject.iCortex).FileName, ...
                      sSubject.Anatomy(sSubject.iAnatomy).FileName);
end

bst_progress('inc',1);
bst_progress('text', 'Saving results');

add_vol_data(sMRI,voronoi, voronoi_fn, ...
             ['Voronoi interpolator for ' sSubject.Anatomy(sSubject.iAnatomy).Comment ...
              ' onto ' sSubject.Surface(sSubject.iCortex).Comment], iSubject);
OutputFiles = {'import'};

bst_progress('stop');

end

function [vol_voro, sMri] = Compute(cortex_file, anatomy_file, segmentation_file)

% Obtain the cortical surface
sCortex = in_tess_bst(cortex_file);

% Obtain the anatomical MRI
sMri = in_mri_bst(anatomy_file);

% Obtain the binary mask around the cortical surface
if isfield(sCortex,'tess2mri_interp') && ~isempty(sCortex.tess2mri_interp)
    tess2mri_interp     = sCortex.tess2mri_interp;
else
    tess2mri_interp     = tess_interp_mri(cortex_file, sMri);
end

index_binary_mask   = (sum(tess2mri_interp,2) >0);
binary_volume       = zeros(size(sMri.Cube));
binary_volume(index_binary_mask) = 1;

% Dilate the mask of the WG interface surface to obtain the gray matter
% (corrected by alexis)
nb_dilations = 5;
binary_volume_dilated = mri_dilate(binary_volume, nb_dilations);

% Convert the coordinates of the cortical surface
Vertices = sCortex.Vertices;

% Convert coordinates
% Vertices: SCS->MRI
Vertices = cs_convert(sMri, 'scs', 'mri', Vertices) * 1000;

% Vertices: MRI(MM)->MRI(Voxels)
Vertices    = bst_bsxfun(@rdivide, Vertices, sMri.Voxsize);
A_nodes     = Vertices'  ;

% Regeneration of a volume with only the seed points 
vol_seeds   = zeros( size(sMri.Cube)  ) ;

% Position the seeds of the Mesh in the MRI voxel grid  : N
N   = zeros(3,size(A_nodes,2));
N1  = zeros(3,size(A_nodes,2));

for i=1:size(A_nodes,2)
    N( :, i ) = A_nodes( :, i ) ; %inv( mat( 1:3, 1:3 ) ) * A_nodes( :, i ) ;
    
    N1(1,i) = round(N(1,i));% (Xmax-N( 1, i )) / vox_size(1) );
    N1(2,i) = round(N(2,i));%round( (Ymax-N( 2, i )) / vox_size(2) );
    N1(3,i) = round(N(3,i));%round( (Zmax-N( 3, i )) / vox_size(3) );
    
    if vol_seeds( N1(1,i), N1( 2,i) , N1(3,i) ) == 0
        vol_seeds(N1(1,i), N1(2,i), N1(3,i)) = i  ;
    end
end

% Compute Voronoi Diagram for interpolation 
xt          = N1';
ListRes     = floor(xt+0.5);
distance    = 'd8';
vox_size    = sMri.Voxsize ;

vol_voro = dg_voronoi(binary_volume_dilated, vox_size, ListRes, distance);

if nargin > 2
    GM_mask = get_grey_matter_mask(segmentation_file);

    tmp = -1 * ones(size(vol_voro));
    tmp(GM_mask) = vol_voro(GM_mask);

    vol_voro = tmp; 

end
end


function GM_mask = get_grey_matter_mask(segmentation_file)

    sSegmentation = in_mri_bst(segmentation_file);
    GM_mask       = false(size(sSegmentation.Cube)); 

    if ~isfield(sSegmentation,'Labels') ||  isempty(sSegmentation.Labels)
        bst_error('Unrecognized segmentation.');
        return;
    end

    idx = [];
    if any(contains({sSegmentation.Labels{:,2}}, 'Cortex'))
        idx = cell2mat({sSegmentation.Labels{contains({sSegmentation.Labels{:,2}}, {'Cortex', 'Hippocampus', 'Thalamus', 'Caudate', 'Putamen', 'Pallidum', 'Amygdala', 'Cerebellum R', 'Cerebellum L'}),1}});
    elseif any(contains({sSegmentation.Labels{:,2}}, 'Gray'))
        idx = cell2mat({sSegmentation.Labels{contains({sSegmentation.Labels{:,2}}, 'Gray'),1}});
    end

    if isempty(idx)
        bst_error('Unrecognized segmentation.');
        return;
    end    

    for i = 1:length(idx)               
        GM_mask(sSegmentation.Cube == idx(i)) = true;
    end

end

function voronoi_fn = get_voronoi_fn(sSubject)
    i_voronoi = find(strcmp( {sSubject.Anatomy.Comment}, sprintf('Voronoi interpolator for %s onto %s',sSubject.Anatomy(sSubject.iAnatomy).Comment, sSubject.Surface(sSubject.iCortex).Comment)));
    
    if any(i_voronoi)
        voronoi_fn = file_fullpath(sSubject.Anatomy(i_voronoi(1)).FileName);
    else
        [cortex_root, cortex_bfn, cortex_ext] = fileparts(sSubject.Surface(sSubject.iCortex).FileName);
        [anat_root, anat_bfn, anat_ext] = fileparts(sSubject.Anatomy(sSubject.iAnatomy).FileName);
        voronoi_bfn = [anat_bfn '_' cortex_bfn  '_voronoi' anat_ext];
        subjectSubDir = bst_fileparts(sSubject.FileName);
        ProtocolInfo = bst_get('ProtocolInfo');
        voronoi_fn =  bst_fullfile(ProtocolInfo.SUBJECTS, subjectSubDir, voronoi_bfn);
    end    
end

function median_volume = get_median_voronoi_volume(sVoronoi)
%==========================================================================
% Returns the median volume of the voronoi cell in mm^3
% Counts how many voxels are in each voronoi cells
%==========================================================================

    counts = groupcounts(sVoronoi.Cube(~isnan(sVoronoi.Cube) & sVoronoi.Cube > 0));    % Volume of one voxel in mm^3
    volume_voxel  = sVoronoi.Voxsize(1) *  sVoronoi.Voxsize(2)  *  sVoronoi.Voxsize(3);    % Median volume of the voronoi cells in mm^3
    median_volume = median(counts) * volume_voxel;
end



function sSubject = add_vol_data(sMri, data, vol_fn, vol_comment, iSubject, history_comment)
%% Save a volume map to brainstorm with given data

ProtocolSubjects = bst_get('ProtocolSubjects');

if (iSubject == 0) % Default subject
	sSubject = ProtocolSubjects.DefaultSubject;
else % Normal subject
    sSubject = ProtocolSubjects.Subject(iSubject);
    if sSubject.UseDefaultAnat
        sSubject = ProtocolSubjects.DefaultSubject;
        iSubject = 0;
    end 
end

if ~exist(vol_fn,'file')
    is_new = 1;
else
    is_new = 0;
end

sMri.Cube       = data;
sMri.Comment    = vol_comment;
sMri.Histogram  = mri_histogram(sMri.Cube);

if nargin > 5
    sMri = bst_history('add', sMri, 'import', history_comment);
end

% Save new MRI in Brainstorm format
out_mri_bst(sMri, vol_fn,'v6');

% If the MRI already exist, no need to go further.
if ~is_new 
    return;
end
%% ===== STORE NEW MRI IN DATABASE ======
% New anatomy structure
iAnatomy = find(cell2mat(cellfun(@(a) ~isempty(a), ...
                   strfind({sSubject.Anatomy.FileName}, vol_fn),...
                   'uni', false )));
% Add new anatomy
if isempty(iAnatomy)
    iAnatomy = length(sSubject.Anatomy) + 1;
end
sSubject.Anatomy(iAnatomy) = db_template('Anatomy');
sSubject.Anatomy(iAnatomy).FileName = file_short(vol_fn);
sSubject.Anatomy(iAnatomy).Comment  = vol_comment;

% == Update database ==
% Default subject
if (iSubject == 0)
	ProtocolSubjects.DefaultSubject = sSubject;
% Normal subject 
else
    ProtocolSubjects.Subject(iSubject) = sSubject;
end
bst_set('ProtocolSubjects', ProtocolSubjects);

%% ===== UPDATE GUI =====
% Refresh tree
panel_protocols('UpdateNode', 'Subject', iSubject);
panel_protocols('SelectNode', [], 'subject', iSubject, -1 );
% Save database
db_save();
end
