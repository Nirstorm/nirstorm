classdef OptimalMontageTest < matlab.unittest.TestCase
    
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
        
        function test_OM_from_cortex(testCase)
            global GlobalData;
            subject_name = bst_create_test_subject();
            
            stg_vertex_id = 19552;
            roi_scout_selection = bst_create_scout(subject_name, 'cortex', 'roi_temporal', stg_vertex_id, 2, 'User scouts');

            extent_cm = 3; % centimeter
            head_vertices = process_nst_cpt_fluences_from_cortex('proj_cortex_scout_to_scalp', ...
                                                                 roi_scout_selection, extent_cm * 0.01);
            wavelengths = 685;
            fluence_dir = cpt_spherical_fluences(roi_scout_selection.sSubject, head_vertices, wavelengths);
            bst_process('CallProcess', ...
                'process_nst_OM_from_cortex', [], [], ...
                'scout_sel_roi', roi_scout_selection, ...
                'cortex_to_scalp_extent', extent_cm, ...
                'condition_name', 'OM_test', ...
                'wavelengths', strjoin(arrayfun(@num2str, wavelengths, 'UniformOutput', false), ','), ...
                'data_source', fluence_dir, ...
                'nb_sources', 1, ...
                'nb_detectors', 2, ...
                'nAdjacentDet', 0, ...
                'exist_weight', 0, ...
                'sep_optode', {[0 55],'mm', 0});
            
            testCase.assertEmpty(GlobalData.lastestFullErrMsg);
            % TODO add checks on nb optodes + separations
        end

        
        function test_OM_from_head(testCase)
            global GlobalData;
            subject_name = bst_create_test_subject();
            head_vertex_idx = 6465; %ASSUME: valid for Colin27 4NIRS default head mesh
            head_scout_selection = bst_create_scout(subject_name, 'scalp', 'OM_temporal', head_vertex_idx, 5, 'User scouts');
            
            stg_vertex_id = 19552; %ASSUME: valid for Colin27 4NIRS default mid mesh
            roi_scout_selection = bst_create_scout(subject_name, 'cortex', 'roi_temporal', stg_vertex_id, 2, 'User scouts');
            
            wavelengths = 685;
            fluence_dir = cpt_spherical_fluences(head_scout_selection.sSubject, head_scout_selection.sScout.Vertices, wavelengths);
            
            bst_process('CallProcess', ...
                'process_nst_OM_from_head', [], [], ...
                'scout_sel_head', head_scout_selection, ...
                'scout_sel_roi', roi_scout_selection, ...
                'condition_name', 'OM_test', ...
                'wavelengths', strjoin(arrayfun(@num2str, wavelengths, 'UniformOutput', false), ','), ...
                'data_source', fluence_dir, ... %TODO
                'nb_sources', 1, ...
                'nb_detectors', 2, ...
                'nAdjacentDet', 0, ...
                'exist_weight', 0, ...
                'sep_optode', {[0 55],'mm', 0});
            testCase.assertEmpty(GlobalData.lastestFullErrMsg);
            % TODO add checks on nb optodes + separations
        end

    end
end

function [subject_name, sSubject] = bst_create_test_subject()
%% Ensure that nst_utest protocol exists
ProtocolName = 'nst_utest';
% Delete existing protocol
db_dir = bst_get('BrainstormDbDir');
gui_brainstorm('DeleteProtocol', ProtocolName);

nst_protocol_dir = fullfile(db_dir, ProtocolName);
if exist(nst_protocol_dir, 'dir')
    rmdir(nst_protocol_dir, 's');
end

% Create new protocol
gui_brainstorm('CreateProtocol', ProtocolName, 0, 0);

subject_name = 'test_subject';
[sSubject, iSubject] = db_add_subject(subject_name);

sTemplates = bst_get('AnatomyDefaults');
iTemplate = strcmpi('Colin27_4NIRS', {sTemplates.Name});
db_set_template(iSubject, sTemplates(iTemplate), 0);
db_save();
sSubject = bst_get('Subject', subject_name);
end

function new_scout_sel = bst_create_scout(subject_name, surface_type, scout_name, scout_seed, scout_size, atlas_name)
assert(isvarname(scout_name));

if strcmp(surface_type, 'scalp')
    s_isurf = 'iScalp';
elseif strcmp(surface_type, 'cortex')
    s_isurf = 'iCortex';
else
    disp('ERROR: surf_type must either be "scalp" or "cortex"');
    bst_error('surf_type must either be "scalp" or "cortex"');
end

if nargin < 6
    atlas_name = 'User scouts';
end

[sSubject, iSubject] = bst_get('Subject', subject_name);

sSurf = in_tess_bst(sSubject.Surface(sSubject.(s_isurf)).FileName);
iatlas = strcmp(atlas_name, {sSurf.Atlas.Name});
if isempty(iatlas)
    error(['Atlas not found:' atlas_name]);
end

assert(scout_seed <= size(sSurf.Vertices, 1));
scout_vertices = [scout_seed];
for isize=1:scout_size
    cur_vertices = scout_vertices;
    for ivertex=1:length(cur_vertices)
        vidx = cur_vertices(ivertex);
        scout_vertices = unique([scout_vertices find(sSurf.VertConn(vidx, :))]);
    end
end
scout_idx = size(sSurf.Atlas(iatlas).Scouts,2) + 1;
new_scout = db_template('Scout');
new_scout.Vertices = scout_vertices;
new_scout.Seed = scout_seed;
new_scout.Color = [0,0,0];
new_scout.Label = scout_name;
sSurf.Atlas(iatlas).Scouts(scout_idx) = new_scout;
bst_save(file_fullpath(sSubject.Surface(sSubject.(s_isurf)).FileName), sSurf, 'v7');
db_save();

new_scout_sel.sSubject = bst_get('Subject', subject_name);
new_scout_sel.isurface = sSubject.(s_isurf);
new_scout_sel.iatlas = iatlas;
new_scout_sel.iscout = scout_idx;
new_scout_sel.sScout = new_scout;
end


function fluence_dir = cpt_spherical_fluences(sSubject, head_scout_vertices, wavelengths, dcut)

fluence_dir = fullfile(nst_get_local_user_dir(), 'fluence', 'MRI__Colin27_4NIRS', 'spherical');
if ~exist(fluence_dir, 'dir')
    mkdir(fluence_dir);
end

if nargin < 4
    dcut = 40; %mm
end

disp(['Start generating missing spherical fluences for ' num2str(length(head_scout_vertices)) ' vertices...']);
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

% Vertices: SCS->MRI, MRI(mm)
head_vertices_mri = cs_convert(sMri, 'scs', 'mri', sHead.Vertices) * 1000; % m to mm

[dimx, dimy, dimz] = size(sMri.Cube);
[xx, yy, zz] = meshgrid(1:dimx, 1:dimy, 1:dimz);
all_coords = [xx(:) yy(:) zz(:)] .* repmat(sMri.Voxsize, dimx*dimy*dimz, 1);
fluence_vol = zeros(dimy, dimx, dimz); %TODO: weird order ... to check

% Save MRI to nifti:
% sVol = sMri;
% out_fn =  fullfile(fluence_dir, 't1.nii');
% out_mri_nii(sMri, out_fn, 'float32');

for ivertex=1:length(head_scout_vertices)
    vertex_id = head_scout_vertices(ivertex);
    fluence_wl1_fn = fullfile(fluence_dir, process_nst_import_head_model('get_fluence_fn', vertex_id, wavelengths(1)));
    if ~exist(fluence_wl1_fn, 'file')
        vertex = head_vertices_mri(vertex_id, :);
        rv = repmat(vertex, size(all_coords, 1), 1);
        cp = rv - all_coords;
        fluence_flat = max(0, 1 - sqrt(sum(cp .* cp, 2)) / dcut);
        fluence_vol(:) = fluence_flat;
        to_save = permute(fluence_vol, [2 3 1]);
        if 0
            vertex_vox = round(vertex ./ sMri.Voxsize);
            figure(); 
            subplot(2,2,1); hold on; 
            imagesc(squeeze(sMri.Cube(:,:,vertex_vox(3)))); 
            ii = imagesc(squeeze(fluence_vol(:,:,vertex_vox(3))));
            set(ii, 'AlphaData', 0.75);
            subplot(2,2,2); imagesc(squeeze(fluence_vol(:,vertex_vox(1),:)));
            subplot(2,2,3); imagesc(squeeze(fluence_vol(vertex_vox(2),:,:)));
            sVol.Comment = '';
            sVol.Cube = to_save;
            sVol.Histogram = [];
            [rr, bfn, ext] = fileparts(fluence_wl1_fn);
            out_fn =  fullfile(fluence_dir, [bfn '.nii']);
            out_mri_nii(sVol, out_fn, 'float32');    
        end
        fluence_flat_sparse_vol = sparse(double(to_save(:)));
        reference_voxel_index = round(vertex ./ sMri.Voxsize);
        save(fluence_wl1_fn, 'fluence_flat_sparse_vol','reference_voxel_index');
        for iwl=2:length(wavelengths)
            wl = wavelengths(iwl);
            fluence_fn = fullfile(fluence_dir, process_nst_import_head_model('get_fluence_fn', vertex_id, wl));
            if ~exist(fluence_fn, 'file')
                copyfile(fluence_wl1_fn, fluence_fn);
            end
        end
    end
end
disp('Done generating missing spherical fluences.');
end

function fluence = load_fluence_spherical(vertex_id, sInputs)

ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
wavelengths = ChannelMat.Nirs.Wavelengths;

fluence_bfn =  sprintf('fluence_spherical_%d.mat', vertex_id); 
fluence_fn = bst_fullfile(bst_get('BrainstormUserDir'), 'defaults', ...
                              'nirstorm', fluence_bfn);
if ~file_exist(fluence_fn)
    % Create folder
    if ~file_exist(bst_fileparts(fluence_fn))
        mkdir(bst_fileparts(fluence_fn));
    end
    
    [sSubject, iSubject] = bst_get('Subject', sInputs.SubjectName);
    anatomy_file = sSubject.Anatomy(sSubject.iAnatomy).FileName;    
    head_mesh_fn = sSubject.Surface(sSubject.iScalp).FileName;
    
    % Obtain the head mesh
    sHead = in_tess_bst(head_mesh_fn);
    
    % Obtain the anatomical MRI
    sMri = in_mri_bst(anatomy_file);
    [dimx, dimy, dimz] = size(sMri.Cube);
    
    vertex = sHead.Vertices(vertex_id, :);
    % Vertices: SCS->MRI, MRI(MM)->MRI(Voxels)
    vertex = round(cs_convert(sMri, 'scs', 'mri', vertex) * 1000 ./ sMri.Voxsize);
    
    dcut = 30; %mm
    svox = sMri.Voxsize;
    fluence_vol = zeros(dimx, dimy, dimz);
    for i=1:dimx
        for j=1:dimy
            for k=1:dimz
                cp = (vertex - [i,j,k]) .* svox;
                fluence_vol(i,j,k) = max(0, 1 - sqrt(sum(cp .* cp)) / dcut);
            end
        end
    end
    save(fluence_fn, 'fluence_vol');
else
    load(fluence_fn, 'fluence_vol');
end

for iwl=1:length(wavelengths)
    fluence{iwl}.fluence.data = fluence_vol;
end
    
end

