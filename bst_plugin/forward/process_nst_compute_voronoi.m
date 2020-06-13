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

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute Voronoi volume-to-cortex interpolator';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS';
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
    
    sProcess.options.segmentation_label.Type    = 'radio_line';
    sProcess.options.segmentation_label.Comment = {'1:skin, 2:skull, 3:CSF, 4:GM, 5:WM', '5: skin,  4: skull, 3: CSF, 2: GM, 1: WM','Segmentation label: '};
    sProcess.options.segmentation_label.Value   = 1;  

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
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
voronoi_fn = get_voronoi_fn(sSubject);

seg_label = 'segmentation_5tissues';
segmentation_id = find(strcmp(seg_label, {sSubject.Anatomy.Comment}));
if ~isempty(segmentation_id)               
    voronoi = Compute(sSubject.Surface(sSubject.iCortex).FileName, ...
                      sSubject.Anatomy(sSubject.iAnatomy).FileName, ...
                      sSubject.Anatomy(segmentation_id).FileName, ...
                      sProcess.options.segmentation_label.Value);
else
    msg = ['MRI segmentation (' seg_label ') not found. ' ...
           'Interpolator cannot be constrained to grey matter, so expect PVE.'];
    disp(['BST Warning> ' msg]);
    bst_report('Warning', sProcess, sInputs, msg);
    voronoi = Compute(sSubject.Surface(sSubject.iCortex).FileName, ...
                      sSubject.Anatomy(sSubject.iAnatomy).FileName,...
                      sProcess.options.segmentation_label.Value);
    if isempty(voronoi)
       return;
    end
end

add_vol_data(voronoi, voronoi_fn, ...
             ['Voronoi interpolator for ' sSubject.Anatomy(sSubject.iAnatomy).Comment ...
              ' onto ' sSubject.Surface(sSubject.iCortex).Comment], iSubject);
OutputFiles = {'import'};
end

function vol_voro = Compute(cortex_file, anatomy_file, segmentation_file,segmentation_label)

if nargin < 4
    segmentation_label=1;
end    
    

% Obtain the cortical surface
sCortex = in_tess_bst(cortex_file);

% Obtain the anatomical MRI
sMri = in_mri_bst(anatomy_file);
[dimx, dimy, dimz] = size(sMri.Cube);


% Obtain the binary mask around the cortical surface
tess2mri_interp = tess_interp_mri(cortex_file, sMri);
index_binary_mask = (sum(tess2mri_interp,2) >0);
binary_volume = zeros(dimx, dimy, dimz);
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
Vertices = bst_bsxfun(@rdivide, Vertices, sMri.Voxsize);

% Generate the seed volume 
[ Xmax, Ymax, Zmax ] = size(sMri.Cube);

A_nodes = [Vertices' ] ;

N1 = [];

% Regeneration of a volume with only the seed points 
vol_seeds = zeros( size(sMri.Cube)  ) ;
Duplicates = {};

% Position the seeds of the Mesh in the MRI voxel grid  : N
for i=1:size(A_nodes,2)
    N( :, i ) = A_nodes( :, i ) ; %inv( mat( 1:3, 1:3 ) ) * A_nodes( :, i ) ;
    
    N1(1,i) = round(N(1,i));% (Xmax-N( 1, i )) / vox_size(1) );
    N1(2,i) = round(N(2,i));%round( (Ymax-N( 2, i )) / vox_size(2) );
    N1(3,i) = round(N(3,i));%round( (Zmax-N( 3, i )) / vox_size(3) );
    
    if vol_seeds( N1(1,i), N1( 2,i) , N1(3,i) ) ~= 0
        Duplicates{end+1} = {vol_seeds( N1(1,i), N1(2,i) , N1(3,i) ), i};
    else
        vol_seeds(N1(1,i), N1(2,i), N1(3,i)) = i  ;
    end
end
vol_seeds( vol_seeds==0 ) = -1 ;

% Compute Voronoi Diagram for interpolation 
xt=N1';
List = floor(xt+0.5);
ListRes = List;
distance='d8';
vox_size = sMri.Voxsize ;
bst_progress('start', 'MRI/Surface Voronoi interpolator','Computing Voronoi partitioning (~3min) ...', 1, 2);
vol_voro = dg_voronoi(binary_volume_dilated, vox_size, ListRes, distance);
% vol_voro = binary_volume_dilated; %HACK
if nargin > 2
    sSegmentation = in_mri_bst(segmentation_file);
    if segmentation_label == 1
        sSegmentation.Cube = nst_prepare_segmentation(sSegmentation.Cube,{1,2,3,4,5});
    elseif segmentation_label == 2
        sSegmentation.Cube = nst_prepare_segmentation(sSegmentation.Cube,{5,4,3,2,1});
    end    
    
    if all(sSegmentation.Histogram.fncX == 0:5)
        vol_voro(sSegmentation.Cube ~= 4) = -1;
    elseif any(sSegmentation.Histogram.fncX == 204)
        vol_voro(sSegmentation.Cube ~= 204) = -1;
    else
        bst_error('Unrecognized segmentation labels: should be from 0 to 5');
        vol_voro = [];
        return;
    end
end

bst_progress('inc',1);
bst_progress('stop');
end


function voronoi_fn = get_voronoi_fn(sSubject)
[cortex_root, cortex_bfn, cortex_ext] = fileparts(sSubject.Surface(sSubject.iCortex).FileName);
[anat_root, anat_bfn, anat_ext] = fileparts(sSubject.Anatomy(sSubject.iAnatomy).FileName);
voronoi_bfn = [anat_bfn '_' cortex_bfn  '_voronoi' anat_ext];
subjectSubDir = bst_fileparts(sSubject.FileName);
ProtocolInfo = bst_get('ProtocolInfo');
voronoi_fn =  bst_fullfile(ProtocolInfo.SUBJECTS, subjectSubDir, voronoi_bfn);
end


function sSubject = add_vol_data(data, vol_fn, vol_comment, iSubject, history_comment)
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

sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
assert(all(size(sMri.Cube) == size(data)));

sMri.Cube = data;
sMri = rmfield(sMri, 'Histogram');
sMri.Comment = vol_comment;

if nargin > 4
    sMri = bst_history('add', sMri, 'import', history_comment);
end

% Save new MRI in Brainstorm format
sMri_out = out_mri_bst(sMri, vol_fn);

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
