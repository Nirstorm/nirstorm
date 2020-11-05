function varargout = process_nst_OM_from_head( varargin )

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
% Authors: Alexis Machado, Thomas Vincent, ZhengChen Cai (2017-2018)

% TODO: expose projection distance from cortex ROI to head
% TODO: only one wavelength (quickly test if montage change from one wl to
%      another)
% TODO: test limitation of model size -> check with 32x32 or 16x16.
% TODO: document computation time
% TODO: doc params: 
%   - nb sources/det is the nb of optodes placed by algo then post proc to
%     remove positions
%   - nb adjacent is MIN (can be higher). Has to be associated to dist MAX
% TODO: fix paper reference + ref to equations

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute optimal montage from head scout';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1408;
sProcess.Description = '';
sProcess.isSeparator = 0;
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

sProcess.options = nst_add_scout_sel_options(struct(), 'head', 'Head scout (search space):', ...
                                             'scalp', {'User scouts'}, 0); 
sProcess.options = add_OM_options(sProcess.options);
end

function options = add_OM_options(options)

% Add selector of cortical scout for target VOI
options = nst_add_scout_sel_options(options, 'roi', 'Cortical scout (target ROI):', ...
                                    'cortex', {'User scouts'}, 1);

options.condition_name.Comment = 'Output condition name:';
options.condition_name.Type = 'text';
options.condition_name.Value = '';

options.segmentation_label.Type    = 'radio_line';
options.segmentation_label.Comment   = {'1:skin, 2:skull, 3:CSF, 4:GM, 5:WM', '5: skin,  4: skull, 3: CSF, 2: GM, 1: WM','Segmentation label: '};
options.segmentation_label.Value   = 1;

options.wavelengths.Comment = 'Wavelengths (nm) [coma-separated list]';
options.wavelengths.Type    = 'text';
options.wavelengths.Value = '';

options.data_source.Comment = 'Fluence Data Source (URL or path)';
options.data_source.Type    = 'text';
options.data_source.Value = [nst_get_repository_url() '/fluence/'];

options.nb_sources.Comment = 'Number of sources:';
options.nb_sources.Type = 'value';
options.nb_sources.Value = {4,'',0};

options.nb_detectors.Comment = 'Number of detectors:';
options.nb_detectors.Type = 'value';
options.nb_detectors.Value =  {8,'',0};

options.nAdjacentDet.Comment = 'Number of Adjacent:';
options.nAdjacentDet.Type = 'value';
options.nAdjacentDet.Value = {2,'',0};

options.sep_optode.Comment = 'Range of optodes distance:';
options.sep_optode.Type = 'range';
options.sep_optode.Value = {[15 55],'mm',0};

options.sepmin_SD.Comment = 'Minimum source detector distance:';
options.sepmin_SD.Type = 'value';
options.sepmin_SD.Value = {15,'mm',0};

options.exist_weight.Comment = 'Use existing weight tables (speed up)';
options.exist_weight.Type = 'checkbox';
options.exist_weight.Value = 0;

SelectOptions = {...
    '', ...                            % Filename
    '', ...                            % FileFormat
    'save', ...                        % Dialog type: {open,save}
    'Select output folder...', ...     % Window title
    'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
    'single', ...                      % Selection mode: {single,multiple}
    'dirs', ...                        % Selection mode: {files,dirs,files_and_dirs}
    {{'.folder'}, '*.*'}, ... % Available file formats
    'MriOut'};                         % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
% Option definition
% TODO: add flag to enable ouput
options.outputdir.Comment = 'Folder for weight table:';
options.outputdir.Type    = 'filename';
options.outputdir.Value   = SelectOptions;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% Get scout vertices & load head mesh
head_scout_selection = nst_get_option_selected_scout(sProcess.options, 'head');
head_vertices = head_scout_selection.sScout.Vertices;

sSubject = head_scout_selection.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);

OutputFiles = Compute(sProcess, sSubject, sHead, head_vertices, sInputs);

end

function OutputFiles = Compute(sProcess, sSubject, sHead, head_vertex_ids, sInputs)
OutputFiles = {};

cplex_url = 'https://www.ibm.com/us-en/marketplace/ibm-ilog-cplex/resources';

try
    cplx = Cplex();
    cplex_version = strsplit(cplx.getVersion(), '.');
    if str2double(cplex_version(1)) < 12 || ...
            (length(cplex_version) > 1 && str2double(cplex_version(1)) < 3)
        bst_error(['CPLEX >12.3 required. See ' cplex_url]);
        return
    end
catch
    bst_error(['CPLEX >12.3 required. See ' cplex_url]);
    return
end

condition_name = sProcess.options.condition_name.Value;
if isempty(condition_name)
    condition_name = 'planning_optimal_montage';
end

% Extract head point coordinates
head_vertices_coords = sHead.Vertices(head_vertex_ids, :);

try
    scan_res = textscan(sProcess.options.wavelengths.Value, '%d,');
catch
    bst_report('Error', sProcess, [], 'List of wavelengths must be integers separated by comas');
    return
end
wavelengths = double(scan_res{1}');

if(length(wavelengths) < 1) % TODO: at least 2 wavelengths ?
    bst_report('Error', sProcess, [], 'List of wavelengths must not be empty');
    return
end

% Obtain the anatomical MRI
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
cubeSize = size(sMri.Cube);

options.nb_sources = sProcess.options.nb_sources.Value{1} ;
options.nb_detectors = sProcess.options.nb_detectors.Value{1};
options.nAdjacentDet = sProcess.options.nAdjacentDet.Value{1};
options.sep_optode_min = sProcess.options.sep_optode.Value{1}(1);
options.sep_optode_max  = sProcess.options.sep_optode.Value{1}(2);
options.sep_SD_min = sProcess.options.sepmin_SD.Value{1};
options.cubeSize = cubeSize;
options.outputdir = sProcess.options.outputdir.Value{1};
options.exist_weight = sProcess.options.exist_weight.Value;

if options.sep_optode_min > options.sep_SD_min
    bst_error(sprintf('ERROR: The minimum distance between source and detector has to be larger than the minimum optodes distance'));
end
% TODO: enable selection of multiple ROIs
% TODO: assert that selected subjecct is the same as for the head scout
roi_scout_selection = nst_get_option_selected_scout(sProcess.options, 'roi');
sCortex = in_tess_bst(roi_scout_selection.sSubject.Surface(roi_scout_selection.isurface).FileName);
sScoutsFinal = {roi_scout_selection.sScout};

vois = cell(length(sScoutsFinal), 1);
for iroi=1:length(sScoutsFinal)
    sROI = extract_scout_surface(sCortex, sScoutsFinal{iroi});
    sROI.Vertices = cs_convert(sMri, 'scs', 'voxel', sROI.Vertices);
    
    % Get VOI from cortical ROI -> interpolate
    tess2mri_interp = tess_tri_interp(sROI.Vertices, sROI.Faces, cubeSize);
    index_binary_mask = (sum(tess2mri_interp,2) >0);
    voi = zeros(cubeSize(1), cubeSize(2), cubeSize(3));
    if 1
        voronoi_fn = process_nst_compute_voronoi('get_voronoi_fn', sSubject);
        
        if ~exist(voronoi_fn, 'file')
            process_nst_compute_voronoi('Run', sProcess, sInputs);%TODO: test this because compute voro expect subject name as input option
        end
        voronoi_bst = in_mri_bst(voronoi_fn);
        voronoi = voronoi_bst.Cube;
        % Load segmentation
        segmentation_name = 'segmentation_5tissues';
        iseg = 0;
        for ianat = 1:size(sSubject.Anatomy,2)
            if any(strcmp(sSubject.Anatomy(ianat).Comment, segmentation_name))
                iseg = ianat;
            end
        end
        if iseg == 0
            bst_error(sprintf('ERROR: Please import segmentation file as MRI and rename it as "%s"', segmentation_name));
        end
        seg = in_mri_bst(sSubject.Anatomy(iseg).FileName);
        if sProcess.options.segmentation_label.Value == 1
            seg.Cube = nst_prepare_segmentation(seg.Cube,{1,2,3,4,5});
        elseif sProcess.options.segmentation_label.Value == 2
            seg.Cube = nst_prepare_segmentation(seg.Cube,{5,4,3,2,1});
        end    
        %TODO: make sure segmentation has proper indexing from 0 to 5
        %      Sometimes goes between 0 and 255 because of encoding issues
        voronoi_mask = (voronoi > -1) & ~isnan(voronoi) & (seg.Cube == 4) & ismember(voronoi,sScoutsFinal{iroi}.Vertices);
%         gm_mask = (seg.Cube == 4);
%         voi_mask = ismember(voronoi,sScoutsFinal{iroi}.Vertices);
        
        % Save MRI to nifti:
%         sVol = sMri;
%         sVol.Cube(:) = gm_mask(:) * 1.0;
%         out_fn =  fullfile('/home/tom/tmp/', 'GM.nii');
%         out_mri_nii(sVol, out_fn, 'float32');
% 
%         sVol = sMri;
%         sVol.Cube = seg.Cube;
%         out_fn =  fullfile('/home/tom/tmp/', 'seg.nii');
%         out_mri_nii(sVol, out_fn, 'float32');
%         
%         
%         sVol = sMri;
%         sVol.Cube(:) = voi_mask(:) * 1.0;
%         out_fn =  fullfile('/home/tom/tmp/', 'voi_mask.nii');
%         out_mri_nii(sVol, out_fn, 'float32');
        
        voi(voronoi_mask) = 1;
    end
    %voi(index_binary_mask) = 1;
    voi_flat = voi(:);
    vois{iroi} = sparse(voi_flat > 0);
    if nnz(voi)==0
        bst_error('VOI mask is empty after intersecting with Voronoi and GM masks');
        return;
    end
end
% Get fluences
sparse_threshold = 1e-6;
% TODO: sparsify saved fluences (faster loading, less memory during load)
fluence_data_source = sProcess.options.data_source.Value;
if ~options.exist_weight
    [fluence_volumes,all_reference_voxels_index] = process_nst_import_head_model('request_fluences', head_vertex_ids, ...
        sSubject.Anatomy(sSubject.iAnatomy).Comment, ...
        wavelengths, fluence_data_source, sparse_threshold);
    if isempty(fluence_volumes)
        return;
    end
end

if options.exist_weight
    [montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords, [], wavelengths, vois, options, []);
else
    [montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords, fluence_volumes, wavelengths, vois, options, all_reference_voxels_index);
end
% montage_pairs contains head vertex indexes -> remap to 1-based consecutive indexes
src_indexes = zeros(max(montage_pairs(:, 1)), 1);
det_indexes = zeros(max(montage_pairs(:, 2)), 1);

% TODO: compute head model

% Forge channels from head points pairs & wavelengths
nChannels = size(montage_pairs, 1) * length(wavelengths);
ChannelMat = db_template('channelmat');
ChannelMat.Comment = 'NIRS-BRS channels';
ChannelMat.Channel = repmat(db_template('channeldesc'), [1, nChannels]);
ChannelMat.Nirs.Wavelengths = wavelengths;
iChan = 1;
det_next_idx = 1;
src_next_idx = 1;
for ipair=1:size(montage_pairs, 1)
    ihead_vertex_src = montage_pairs(ipair, 1);
    ihead_vertex_det = montage_pairs(ipair, 2);
    if src_indexes(ihead_vertex_src) == 0
        idx_src = src_next_idx;
        src_indexes(ihead_vertex_src) = src_next_idx;
        src_next_idx = src_next_idx + 1;
    else
        idx_src = src_indexes(ihead_vertex_src);
    end
    if det_indexes(ihead_vertex_det) == 0
        idx_det = det_next_idx;
        det_indexes(ihead_vertex_det) = det_next_idx;
        det_next_idx = det_next_idx + 1;
    else
        idx_det = det_indexes(ihead_vertex_det);
    end
    disp(['Channel S', num2str(idx_src), 'D' num2str(idx_det) ' >>> Distance: ',...
        num2str(round(pdist2(head_vertices_coords(ihead_vertex_src, :),head_vertices_coords(ihead_vertex_det, :)).*1000,1)), 'mm    ', 'Weight: ',...
        num2str(round(montage_weight(ipair,:),3))]);
    for iwl=1:length(wavelengths)
        
        ChannelMat.Channel(iChan).Name    = sprintf('S%dD%dWL%d', idx_src, idx_det, ...
            ChannelMat.Nirs.Wavelengths(iwl));
        ChannelMat.Channel(iChan).Type    = 'NIRS';
        
        ChannelMat.Channel(iChan).Loc(:,1)  = head_vertices_coords(ihead_vertex_src, :);
        ChannelMat.Channel(iChan).Loc(:,2)  = head_vertices_coords(ihead_vertex_det, :);
        ChannelMat.Channel(iChan).Orient  = [];
        ChannelMat.Channel(iChan).Weight  = 1;
        ChannelMat.Channel(iChan).Comment = [];
        ChannelMat.Channel(iChan).Group = sprintf('WL%d', round(ChannelMat.Nirs.Wavelengths(iwl)));
        
        % OM from cap - make the link between optodes and references name
        % in eeg cap 
        if isfield(sProcess,'reference') && isfield(sProcess,'additional_channel')
            idx_eeg_src= knnsearch(sProcess.reference.loc, head_vertices_coords(ihead_vertex_src, :));
            idx_eeg_det= knnsearch(sProcess.reference.loc, head_vertices_coords(ihead_vertex_det, :));
            
            ChannelMat.Channel(iChan).Comment =  sprintf('%s - %s',sProcess.reference.name{idx_eeg_src},sProcess.reference.name{idx_eeg_det} ); 
            sProcess.additional_channel.Channel(idx_eeg_src).Comment=sprintf('S%d', idx_src);
            sProcess.additional_channel.Channel(idx_eeg_det).Comment=sprintf('D%d', idx_det);
        end     
        iChan = iChan + 1;
    end
end
if isfield(sProcess,'reference') &&  isfield(sProcess,'additional_channel')
    for ieeg=1:size(sProcess.additional_channel.Channel,2)
        if ~isempty( sProcess.additional_channel.Channel(ieeg).Comment)
            ChannelMat.Channel(iChan).Name    = sProcess.additional_channel.Channel(ieeg).Name;
            ChannelMat.Channel(iChan).Type    = sProcess.additional_channel.Channel(ieeg).Type;
            ChannelMat.Channel(iChan).Loc     = sProcess.additional_channel.Channel(ieeg).Loc;
            ChannelMat.Channel(iChan).Comment = sProcess.additional_channel.Channel(ieeg).Comment;
            ChannelMat.Channel(iChan).Orient  = [];
            ChannelMat.Channel(iChan).Weight  = 1;
        
            iChan = iChan + 1;
        end
    end 
end

iStudy = db_add_condition(sSubject.Name, condition_name);
sStudy = bst_get('Study', iStudy);
db_set_channel(iStudy, ChannelMat, 1, 0);
end

function [montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords, fluence_volumes, ...
                                                                  wavelengths, vois, options, all_reference_voxels_index)

% TODO: filter holders by distance to VOI
% TODO: subsample holder mesh if head mesh edge are too small
%       -> maybe build head_OM with fewer vertices?

if ~isfield(options, 'nb_sources')
    bst_error('Number of sources is required');
end

if ~isfield(options, 'nb_detectors')
    bst_error('Number of detectors is required');
end

if ~isfield(options, 'nAdjacentDet')
    bst_error('Number of adjacent is required');
end

if ~isfield(options, 'sep_optode_min')
    bst_error('Minimum distance of optodes is required');
end
if ~isfield(options, 'sep_optode_max')
    bst_error('Maximum distance of optodes is required');
end

holder_distances = pdist2(head_vertices_coords, head_vertices_coords).*1000; % mm
nHolders = size(head_vertices_coords, 1);
nVois = length(vois);
iwl = 1;
weight_tables = cell(nVois, 1);
mri_zeros = zeros(options.cubeSize);
if ~options.exist_weight
    bst_progress('start', 'Compute weights','Computing summed sensitivities of holder pairs...', 1, nVois * nHolders^2);
    for ivoi=1:nVois
        weight_tables{ivoi} = sparse(nHolders, nHolders);
        for isrc=1:nHolders
            for idet=1:nHolders
                if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet)< options.sep_optode_max
                    %A=normFactor*fluenceSrc.*fluenceDet./diff_mask(idx_vox);
                    fluenceSrc = full(fluence_volumes{isrc}{iwl}(vois{ivoi}));
                    fluenceDet = full(fluence_volumes{idet}{iwl}(vois{ivoi}));
                    ref_det_pos = sub2ind(options.cubeSize, all_reference_voxels_index{idet}{iwl}(1), all_reference_voxels_index{idet}{iwl}(2), all_reference_voxels_index{idet}{iwl}(3));
                    %fluenceSD = mri_zeros;
                    %fluenceSD(:) = fluence_volumes{isrc}{iwl};
                    if full(fluence_volumes{isrc}{iwl}(ref_det_pos))==0
                        sensitivity = 0;
                    else
                        sensitivity = fluenceSrc .* fluenceDet ./ fluence_volumes{isrc}{iwl}(ref_det_pos); % Asymmetry
                        % bst_progress('text', ['Maximum sensitivity: ' num2str(round(max(sensitivity),5))]);
                    end
                    weight_tables{ivoi}(isrc, idet) = sum(sensitivity(:));
                end
            end
            bst_progress('inc', nHolders);
        end
        if(nnz(weight_tables{ivoi}) == 0)
            bst_error(sprintf('Weight table is null for VOI %d', ivoi));
            return
        end
    end
    bst_progress('stop');
    if ~isempty(options.outputdir)
        save(fullfile(options.outputdir, 'weight_tables.mat'), 'weight_tables');
    end
else
    load (fullfile(options.outputdir, 'weight_tables.mat'),'weight_tables');
end



nS = options.nb_sources; % number of sources
nD = options.nb_detectors; % number of detectors

flag_sep_optode_optode=1;
thresh_sep_optode_optode=[options.sep_optode_min options.sep_optode_max]; % mm (1) no optodes sep below this thresh (2) Optodes above this thresholds do not form pairs

flag_sep_src_det=1;
thresh_min_sep_src_det = options.sep_SD_min; % mm (1) min sep between a source and a det %TODO: why not using thresh_sep_optode_optode(1)?

flag_adjacency=1;
nAdjacentDet= options.nAdjacentDet; % minimal number of adjacent detector in the given range

flag_init=0; % provide an initial solution

%TODO: fix loop over VOIs!

nH = nHolders; % number of holders
xp=zeros(nH,1); nX=size(xp,1); % binary
yq=zeros(nH,1); nY=size(yq,1); % binary
wq_V=zeros(nH,1); nW_V=size(wq_V,1); % float
v=[xp;yq;wq_V]; nVar=size(v,1);

%..........................................................................
% Aeq_1: Equation 11 Sum(xp)=nSrc
%..........................................................................
Aeq_1=zeros(1,nVar);
Aeq_1(1,1:nH)=1;
E_1=nS;

%..........................................................................
% Aeq_2: Equation 11 Sum(yq)=nDet
%..........................................................................
Aeq_2=zeros(1,nVar);
Aeq_2(1,nH+1:2*nH)=1;
E_2=nD;

%..........................................................................
% Display Equality matrix
%..........................................................................
if 0
    matrixes2show={Aeq_1,Aeq_2,[Aeq_1;Aeq_2]};
    for ii=1:2
        figure ('name','equality matrix')
        mat2show=matrixes2show{ii} ;
        colormap(gray(3))
        imagesc([mat2show]);
        %set(gca,'PlotBoxAspectRatio',[size(mat2show,1) size(mat2show,2) 1])
        ylim([0.5 size([mat2show],1)+0.5])
    end
    clear ii
end


 
%==========================================================================
% Create inequality matrixes
%==========================================================================
M=realmax;
%..........................................................................
% Aineq1: Equation 12 : wq_V-Myq<=0
%..........................................................................
Aineq_1=zeros(nH,nVar);
for iH=1:nH
    Aineq_1(iH,nH+iH)=-M;
    Aineq_1(iH,2*nH+iH)=1;
end
I_1=zeros(nH,1);


ipair = 1;
for iVOI=1:nVois
    bst_progress('start', 'Optimization setup', 'Setting up optimization...', 1, 4);
    weight_table = weight_tables{iVOI};
    % eq 27 ??
    
    if 0 weight_table=weight_table/max(weight_table(:)); end % normalize weigth table
    
    %..........................................................................
    % Aineq2: Equation 13: wq_V-Sum(Vpq*xp)<=0
    %..........................................................................
    Aineq_2=zeros(nH,nVar);
    for iH=1:nH
        Aineq_2(iH,2*nH+iH) = 1;
        Aineq_2(iH,1:nH) = -weight_table(:, iH); %Note: this takes into account weight assymetry
    end
    I_2=zeros(nH,1);
    
    %..........................................................................
    % Aineq3: Optode/Optode minimal separation constraints
    %..........................................................................
    if flag_sep_optode_optode
        separations=holder_distances;
        forbid_sets=zeros(nH);
        forbid_sets(separations<thresh_sep_optode_optode(1))=1; %TODO expose as separate parameter
        
        M2=nS+nD;
        Aineq_3=sparse(2*nH,nVar);
        Aineq_3(1:nH,1:2*nH)=[(M2-1)*eye(nH) + forbid_sets  forbid_sets] ;
        Aineq_3(nH +1:2*nH,1:2*nH)=[forbid_sets (M2-1)*eye(nH) + forbid_sets] ;
        I_3=M2*ones(2*nH,1);
        clear M2
        clear forbid_sets separations
    else
        Aineq_3=[];
        I_3=[];
    end
    bst_progress('inc', 1);
    %..........................................................................
    % Aineq4: constraint source/det min separations
    %..........................................................................
    if flag_sep_src_det
        separations=holder_distances;
        forbid_sets=zeros(nH);
        forbid_sets(separations<thresh_min_sep_src_det)=1;
        
        M2=nD;
        Aineq_4=sparse(nH,nVar);
        Aineq_4(1:nH,1:2*nH)=[M2*eye(nH) forbid_sets] ;
        I_4=M2*ones(nH,1);
        clear M2
        clear forbid_sets separations
    else
        Aineq_4=[];
        I_4=[];
    end
    
    %..........................................................................
    % Aineq5: Adjacency constraints
    %..........................................................................
    if flag_adjacency
        separations=holder_distances;
        Allow_sets=zeros(nH);
        Allow_sets(separations>=thresh_sep_optode_optode(1) & separations<=thresh_sep_optode_optode(2))=1;
        
        Aineq_5=sparse(nH,nVar);
        Aineq_5(1:nH,1:2*nH)=[eye(nH)*nAdjacentDet , -Allow_sets];
        I_5=zeros(nH,1);
    else
        Aineq_5=[];
        I_5=[];
    end
    
    %..........................................................................
    % Display InEquality matrix
    %..........................................................................
    if 0
        matrixes2show={Aineq_1,Aineq_2,Aineq_3, Aineq_4,[Aineq_1;Aineq_2;Aineq_3; Aineq_4; Aineq_5]};
        for ii=1:5
            figure('name','inequality matrix')
            %colormap(flipud(autumn(3)))
            mat2show=matrixes2show{ii};
            imagesc(mat2show);
            %set(gca,'PlotBoxAspectRatio',[size(mat2show,1) size(mat2show,2) 1])
            ylim([0.5 size(mat2show,1)+0.5])
        end
        clear ii
    end
    
    bst_progress('inc', 1);
    %==========================================================================
    % initial condition: find sources whose pairs have maximum energy then
    % complete with detectors
    %==========================================================================
    xp_0=zeros(nH,1);
    yq_0=zeros(nH,1);
    [~,idx] = sort(weight_table(:),'descend');
    % Seems to only needed if holder list is not the same size as
    % weight_table
    if 0
        inc=1; src_pos_idx=[];det_pos_idx=[];
        while numel(src_pos_idx)~=nS
            [ia,~]=ind2sub(size(weight_table),idx(inc));
            if ~ismember(ia,src_pos_idx) && ismember(ia,find(ismember(mont.list_hold,holderList)))
                src_pos_idx = [src_pos_idx ia];
                inc=inc+1;
            else
                inc=inc+1;
                continue
            end
        end
        clear inc ia ib
    end
    [src_pos_idx, ~] = ind2sub(size(weight_table), idx(1:nS));
    
    red_weight_table = weight_table(src_pos_idx,:);
    [~,idx] = sort(red_weight_table(:),'descend');
    if 0
        inc=1;
        while numel(det_pos_idx)~=nD
            [~,ib]=ind2sub(size(red_weight_table),idx(inc));
            if   ~ismember(ib,src_pos_idx) && ismember(ib,find(ismember(mont.list_hold,holderList)))
                det_pos_idx=unique([det_pos_idx ib]);
                inc=inc+1;
            else
                inc=inc+1;
                continue
            end
        end
        clear inc ia ib
    end
    bst_progress('inc', 1);
    [~, det_pos_idx] = ind2sub(size(red_weight_table), idx(1:nD));
    
    xp_0(src_pos_idx)=1;
    yq_0(det_pos_idx)=1;
    
    wq_V_O=zeros(nH,1);
    for ii=1:numel(yq_0)
        if yq_0(ii)==1
            wq_V_0(ii) = sum(weight_table(src_pos_idx,ii));
        end
    end
    incumbent_x0=[];
    incumbent_x0=full(sum(reshape(weight_table(src_pos_idx,det_pos_idx),[],1)));
    display(sprintf('incumbent=%g',incumbent_x0))
    clear src_pos_idx det_pos_idx
    x0=[xp_0;yq_0; wq_V_O];
    
    Aeq=[Aeq_1;Aeq_2];
    E=[E_1;E_2];
    
    Aineq=[Aineq_1;Aineq_2;Aineq_3;Aineq_4,;Aineq_5];
    I=[I_1;I_2;I_3;I_4;I_5];
    
    f=[zeros(1,2*nH) ones(1,nH)]';
    lb=[];%[zeros(1,2*nH) zeros(1,nH)]';
    ub=[];%[ones(1,2*nH)  realmax.*ones(1,nH)]';
    ctype=[repmat('B',1,2*nH) repmat('S',1,nH)];
    
    prob = cplexcreateprob('cplexmilp');
    prob.f = f;
    prob.lb = lb; prob.ub = ub;
    prob.ctype=ctype;
    prob.Aineq=Aineq; prob.bineq = I;
    prob.Aeq=Aeq; prob.beq = E;
    prob.x0=[];
    prob.options = [];
    bst_progress('inc', 1);
    bst_progress('stop');
    
    bst_progress('start', 'Optimization','Running optimization with Cplex. May take several minutes (see matlab console) ...', 1, 2);
    
    cplex=Cplex(prob);
    cplex.Model.sense = 'maximize';
    cplex.Param.timelimit.Cur=300;
    
    if flag_init
        cplex.MipStart(1).name='optim';
        cplex.MipStart(1).effortlevel=1;
        cplex.MipStart(1).x=x0;
        cplex.MipStart(1).xindices=int32([1:numel(x0)]');
    end
    results = cplex.solve;
    bst_progress('inc', 2);
    bst_progress('stop');
    
    if ~isfield(results, 'x')
        bst_error(['OM computation failed  at Cplex step:', results.statusstring]);
        return;
    end
    x=results.x;
    x=round(x);
    isources = find(x(1:nH)==1);
    idetectors = find(x(nH+1:2*nH)==1);
    
    for isrc=1:length(isources)
        for idet=1:length(idetectors)
            if holder_distances(isources(isrc),idetectors(idet)) > thresh_sep_optode_optode(1) && ...
                    holder_distances(isources(isrc),idetectors(idet)) < thresh_sep_optode_optode(2) && ...
                    full(weight_table(isources(isrc),idetectors(idet)))
                montage_pairs(ipair,:) = [isources(isrc) idetectors(idet)];
                montage_weight(ipair,:) = full(weight_table(isources(isrc),idetectors(idet)));
                %             disp(['Channel S', num2str(isrc), 'D' num2str(idet) ' >>> Distance: ',...
                %                 num2str(round(holder_distances(isources(isrc),idetectors(idet)),1)), 'mm    ', 'Weight: ',...
                %                 num2str(round(montage_weight(ipair,:),3))]);
                ipair = ipair + 1;
            end
        end
    end
    
    %     sensitivity_thresh = 0.05;
    %
    %montage_pairs = montage_pairs(find(montage_weight > 0),:);
    
    
end

if 0
    %% Dummy positioning
    sources = knnsearch(head_vertices_coords, mean(head_vertices_coords,1));
    nb_dets = 2;
    detectors = randsample(size(head_vertices_coords, 1), nb_dets);
    
    montage_pairs = zeros(length(sources) * length(detectors), 2);
    ipair = 1;
    for isrc=1:length(sources)
        for idet=1:length(detectors)
            montage_pairs(ipair,:) = [sources(isrc) detectors(idet)];
            ipair = ipair + 1;
        end
    end
end
end

function sSurfNew = extract_scout_surface(sSurf, sScouts)

iRemoveVert = setdiff(1:size(sSurf.Vertices,1), [sScouts.Vertices]);
tag = 'OM_scout_extract';

% Unload everything
bst_memory('UnloadAll', 'Forced');

% === REMOVE VERTICES FROM SURFACE FILE ===
% Remove vertices
[Vertices, Faces, Atlas] = tess_remove_vert(sSurf.Vertices, sSurf.Faces, iRemoveVert, sSurf.Atlas);
% Remove the handles of the scouts
for iAtlas = 1:length(Atlas)
    for is = 1:length(Atlas(iAtlas).Scouts)
        Atlas(iAtlas).Scouts(is).Handles = [];
    end
    %         if isfield(Atlas(iAtlas).Scouts, 'Handles');
    %             Atlas(iAtlas).Scouts = rmfield(Atlas(iAtlas).Scouts, 'Handles');
    %         end
end
% Build new surface
sSurfNew = db_template('surfacemat');
sSurfNew.Comment  = [sSurf.Comment '_' tag];
sSurfNew.Vertices = Vertices;
sSurfNew.Faces    = Faces;
sSurfNew.Atlas    = Atlas;
sSurfNew.iAtlas   = sSurf.iAtlas;

% % === SAVE NEW FILE ===
% % Output filename
% NewTessFile = strrep(file_fullpath(sSurf.FileName), '.mat', ['_' tag '.mat']);
% NewTessFile = file_unique(NewTessFile);
% % Save file back
% bst_save(NewTessFile, sSurfNew, 'v7');
% % Get subject
% [sSubject, iSubject] = bst_get('SurfaceFile', sSurf.FileName);
% % Register this file in Brainstorm database
% db_add_surface(iSubject, NewTessFile, sSurfNew.Comment);

end