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
        to_save = permute(fluence_vol, [2 1 3]);
        
        [vmax, imax] = max(to_save(:));
        [vvox_i, vvox_j, vvox_k] = ind2sub([dimx, dimy, dimz], imax);
        vertex_vox = round(vertex ./ sMri.Voxsize);
        assert(all([vvox_i, vvox_j, vvox_k] == vertex_vox));

        if 0 
            figure(); 
            
            subplot(2,2,1); hold on; 
            imagesc(double(squeeze(sMri.Cube(:,:,vertex_vox(3)))) .* (squeeze(to_save(:,:,vertex_vox(3)))*10+1)); 
            
            subplot(2,2,2); imagesc(double(squeeze(sMri.Cube(:,vertex_vox(2),:))) .* (squeeze(to_save(:,vertex_vox(2),:))*10+1));
            subplot(2,2,3); imagesc(double(squeeze(sMri.Cube(vertex_vox(1),:,:))) .* (squeeze(to_save(vertex_vox(1),:,:))*10+1));
            colormap gray;
            
            title(sprintf('vertex %d', vertex_id));
            
            sVol = sMri;
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
    end
    for iwl=2:length(wavelengths)
        wl = wavelengths(iwl);
        fluence_fn = fullfile(fluence_dir, process_nst_import_head_model('get_fluence_fn', vertex_id, wl));
        if ~exist(fluence_fn, 'file')
            copyfile(fluence_wl1_fn, fluence_fn);
        end
    end
end
disp('Done generating missing spherical fluences.');
end
