function varargout = process_nst_OM( varargin )

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
% Authors: 
% Edouard Delaire (2021)
% Thomas Vincent, Alexis Machado, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute optimal montage';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1406;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

% Subject name
sProcess.options.subjectname.Comment = 'Subject name:';
sProcess.options.subjectname.Type    = 'subjectname';
sProcess.options.subjectname.Value   = '';

sProcess.options.fluencesCond.Comment = {'panel_nst_OM', 'Optimal montage parameters'};
sProcess.options.fluencesCond.Type    = 'editpref';
sProcess.options.fluencesCond.Value   = [];
end


function s = str_pad(s,padsize)
    if (length(s) < padsize)
        s = [repmat('&nbsp;', 1, padsize - length(s)), s];
    end
    s = ['<FONT FACE="monospace">' s '</FONT>'];
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFile = {};

if bst_iscompiled()
    bst_error('Optimum montage is not available in the compiled version of brainstorm');
    return;
end        

cplex_url = 'https://www.ibm.com/us-en/marketplace/ibm-ilog-cplex/resources';
try
    cplx = Cplex();
    cplex_version = nst_strsplit(cplx.getVersion(), '.');
    if str2double(cplex_version(1)) < 12 || ...
            (length(cplex_version) > 1 && str2double(cplex_version(1)) < 3)
        bst_error(['CPLEX >12.3 required. See ' cplex_url]);
        return
    end
catch
    bst_error(['CPLEX >12.3 required. See ' cplex_url]);
    return
end

options = sProcess.options.fluencesCond.Value;
    
SubjectName =  options.SubjectName;
sProcess.options.subjectname.Value = SubjectName;
sSubject = bst_get('Subject',SubjectName);

if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
    bst_error('No available Cortex and Head surface for this subject.');
    return;
end    

sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

i_atlas_cortex = strcmp({sCortex.Atlas.Name},options.Atlas_cortex);
i_scout_cortex = strcmp({sCortex.Atlas(i_atlas_cortex).Scouts.Label}, options.ROI_cortex);
ROI_cortex = sCortex.Atlas(i_atlas_cortex).Scouts(i_scout_cortex);   

if isempty(options.Atlas_head) && isempty(options.ROI_head)
    
    cortex_to_scalp_extent = options.Extent;
    cortex_scout.sSubject = sSubject;
    cortex_scout.sScout   = ROI_cortex;

    [head_vertex_ids, ~, ~] = proj_cortex_scout_to_scalp(cortex_scout,cortex_to_scalp_extent.*0.01, 1); 
else    
    i_atlas_head = strcmp({sHead.Atlas.Name},optionsAtlas_head);
    i_scout_head = strcmp({sHead.Atlas(i_atlas_head).Scouts.Label}, optionsROI_head);
    head_vertex_ids = sHead.Atlas(i_atlas_head).Scouts(i_scout_head).Vertices;   
end


condition_name = optionscondition_name;
if isempty(condition_name)
    condition_name = 'planning_optimal_montage';
end

% Extract head point coordinates
head_vertices_coords = sHead.Vertices(head_vertex_ids, :);

try
    wavelengths = str2double(strsplit(optionswavelengths,','));
catch
    bst_report('Error', sProcess, [], 'List of wavelengths must be integers separated by comas');
    return
end

if(length(wavelengths) < 1) % TODO: at least 2 wavelengths ?
    bst_report('Error', sProcess, [], 'List of wavelengths must not be empty');
    return
end

% Obtain the anatomical MRI
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
cubeSize = size(sMri.Cube);

options.nb_sources = optionsnb_sources ;
options.nb_detectors = optionsnb_detectors;
options.nAdjacentDet = optionsnAdjacentDet;
options.sep_optode_min = optionssep_optode(1);
options.sep_optode_max  = optionssep_optode(2);
options.sep_SD_min = optionssepmin_SD;
options.cubeSize = cubeSize;
options.outputdir = optionsoutputdir;
options.exist_weight = optionsexist_weight;

if options.sep_optode_min > options.sep_SD_min
    bst_error(sprintf('ERROR: The minimum distance between source and detector has to be larger than the minimum optodes distance'));
    return;
end


sScoutsFinal = ROI_cortex;

sROI = extract_scout_surface(sCortex, sScoutsFinal);
sROI.Vertices = cs_convert(sMri, 'scs', 'voxel', sROI.Vertices');
    
% Get VOI from cortical ROI -> interpolate
tess2mri_interp = tess_tri_interp(sROI.Vertices, sROI.Faces, cubeSize);
index_binary_mask = (sum(tess2mri_interp,2) >0);

voi = zeros(cubeSize(1), cubeSize(2), cubeSize(3));
voronoi_fn = process_nst_compute_voronoi('get_voronoi_fn', sSubject);

if ~exist(voronoi_fn, 'file')
    bst_process('CallProcess', 'process_nst_compute_voronoi', [], [], ...
                'subjectname',        sSubject.Name, ...
                'segmentation_label', optionssegmentation_label);
end
    
    
voronoi_bst = in_mri_bst(voronoi_fn);
voronoi = voronoi_bst.Cube;
    
% Load segmentation
segmentation_name = 'segmentation_5tissues';
iseg = find(strcmp({sSubject.Anatomy.Comment},segmentation_name));
if isempty(iseg)
    bst_error(sprintf('ERROR: Please import segmentation file as MRI and rename it as "%s"', segmentation_name));
    return;
end
seg = in_mri_bst(sSubject.Anatomy(iseg).FileName);

if  optionssegmentation_label == 1 % GM label = 4
    voronoi_mask = (voronoi > -1) & ~isnan(voronoi) & (seg.Cube == 4) & ismember(voronoi,sScoutsFinal(iroi).Vertices);
elseif  optionssegmentation_label == 2 % GM label = 2
    voronoi_mask = (voronoi > -1) & ~isnan(voronoi) & (seg.Cube == 2) & ismember(voronoi,sScoutsFinal(iroi).Vertices);
end    

voi(voronoi_mask) = 1;
voi_flat = voi(:);
vois = sparse(voi_flat > 0);
if nnz(voi)==0
    bst_error('VOI mask is empty after intersecting with Voronoi and GM masks');
    return;
end

if options.exist_weight && exist(fullfile(options.outputdir , 'weight_tables.mat'))
    tmp=load (fullfile(options.outputdir, 'weight_tables.mat'),'weight_tables');
    options.weight_tables = tmp.weight_tables;
else
    % Compute weight table
    sparse_threshold = 1e-6;
    fluence_data_source = options.data_source;
    [fluence_volumes,all_reference_voxels_index] = process_nst_import_head_model('request_fluences', head_vertex_ids, ...
        sSubject.Anatomy(sSubject.iAnatomy).Comment, ...
        wavelengths, fluence_data_source, sparse_threshold);
    if isempty(fluence_volumes)
        return;
    end
    weight_tables = compute_weights(fluence_volumes,vois, head_vertices_coords,all_reference_voxels_index, options);
    if(nnz(weight_tables) == 0)
        bst_error(sprintf('Weight table is null for VOI %s', sScoutsFinal.label));
        return
    end
    options.weight_tables = weight_tables;
    if ~isempty(options.outputdir)
        save(fullfile(options.outputdir, 'weight_tables.mat'), 'weight_tables');
    end
end

[montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords,options);
% montage_pairs contains head vertex indexes -> remap to 1-based consecutive indexes
src_indexes = zeros(max(montage_pairs(:, 1)), 1);
det_indexes = zeros(max(montage_pairs(:, 2)), 1);


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
        
        iChan = iChan + 1;
    end
end


iStudy = db_add_condition(sSubject.Name, condition_name);
sStudy = bst_get('Study', iStudy);
db_set_channel(iStudy, ChannelMat, 1, 0);

separations = process_nst_separations('Compute',ChannelMat.Channel) * 100; %convert to cm
% Save time-series data
sDataOut              = db_template('data');
sDataOut.F            = separations;
sDataOut.Comment      = 'Separations';
sDataOut.ChannelFlag  = ones(length(separations),1);
sDataOut.Time         = [1];
sDataOut.DataType     = 'recordings';
sDataOut.nAvg         = 1;
sDataOut.DisplayUnits = 'cm';

% Generate a new file name in the same folder
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_chan_dist');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(iStudy, OutputFile, sDataOut);
end

function weight_tables = compute_weights(fluence_volumes,vois, head_vertices_coords,all_reference_voxels_index, options)
    holder_distances = pdist2(head_vertices_coords, head_vertices_coords).*1000; % mm
    nHolders = size(head_vertices_coords, 1);
    iwl = 1;
    weight_tables = sparse(nHolders, nHolders);
    tic
    bst_progress('start', 'Compute weights','Computing summed sensitivities of holder pairs...', 1, nHolders^2);
    for isrc=1:nHolders
        fluenceSrc = full(fluence_volumes{isrc}{iwl}(vois));
        for idet=1:isrc
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet)< options.sep_optode_max
                %A=normFactor*fluenceSrc.*fluenceDet./diff_mask(idx_vox);

                fluenceDet = full(fluence_volumes{idet}{iwl}(vois{ivoi}));
                ref_det_pos = sub2ind(options.cubeSize, all_reference_voxels_index{idet}{iwl}(1), all_reference_voxels_index{idet}{iwl}(2), all_reference_voxels_index{idet}{iwl}(3));
                if full(fluence_volumes{isrc}{iwl}(ref_det_pos))~=0
                    % sensitivity = fluenceSrc .* fluenceDet ./ fluence_volumes{isrc}{iwl}(ref_det_pos); % Asymmetry
                    % ED: fluenceSrc' * fluenceDet is slighty faster than sum ( fluenceSrc .* fluenceDet )
                    sensitivity = fluenceSrc' * fluenceDet; % ./ fluence_volumes{isrc}{iwl}(ref_det_pos); % Asymmetry

                    weight_tables(isrc, idet) = sensitivity; 
                    if isrc ~= idet
                        weight_tables(idet, isrc) = sensitivity; % The matrix is symetrical before the normalization
                    end    
                end
            end    
        end
        bst_progress('inc', nHolders);
    end
    bst_progress('start', 'Compute weights','Normalizing the sensitivities of holder pairs...', 1,  nHolders^2);
    for isrc=1:nHolders
        for idet=1:nHolders
            if holder_distances(isrc, idet) > options.sep_SD_min && holder_distances(isrc, idet)< options.sep_optode_max                   
                ref_det_pos = sub2ind(options.cubeSize, all_reference_voxels_index{idet}{iwl}(1), all_reference_voxels_index{idet}{iwl}(2), all_reference_voxels_index{idet}{iwl}(3));
                if full(fluence_volumes{isrc}{iwl}(ref_det_pos)) ~=0
                    weight_tables(isrc, idet) = weight_tables(isrc, idet) / fluence_volumes{isrc}{iwl}(ref_det_pos); % Asymmetry
                end
            end
        end
        bst_progress('inc', nHolders);
    end
    toc

    bst_progress('stop');
end


function [montage_pairs,montage_weight] = compute_optimal_montage(head_vertices_coords,options)



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

if ~isfield(options, 'weight_tables')
    bst_error('Weight_tables is required');
end

holder_distances = pdist2(head_vertices_coords, head_vertices_coords).*1000; % mm
nHolders = size(head_vertices_coords, 1);

weight_table = options.weight_tables;


nS = options.nb_sources; % number of sources
nD = options.nb_detectors; % number of detectors

flag_sep_optode_optode=1;
thresh_sep_optode_optode=[options.sep_optode_min options.sep_optode_max]; % mm (1) no optodes sep below this thresh (2) Optodes above this thresholds do not form pairs

flag_sep_src_det=1;
thresh_min_sep_src_det = options.sep_SD_min; % mm (1) min sep between a source and a det %TODO: why not using thresh_sep_optode_optode(1)?

flag_adjacency=1;
nAdjacentDet= options.nAdjacentDet; % minimal number of adjacent detector in the given range

flag_init=0; % provide an initial solution

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
bst_progress('start', 'Optimization setup', 'Setting up optimization...', 1, 4);
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
            ipair = ipair + 1;
        end
    end
end

%sensitivity_thresh = 0.05;
%montage_pairs = montage_pairs(find(montage_weight > 0),:);
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


function [head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout, extent_m, save_in_db)

if nargin < 3
    save_in_db = 0;
end
sSubject = cortex_scout.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
dis2head = pdist2(sHead.Vertices, sCortex.Vertices(cortex_scout.sScout.Vertices,:));
head_vertices = find(min(dis2head,[],2) < extent_m); 

% TODO: properly select atlas
exclude_scout = sHead.Atlas.Scouts(strcmp('FluenceExclude', {sHead.Atlas.Scouts.Label}));
if ~isempty(exclude_scout)
    head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
end

limiting_scout = sHead.Atlas.Scouts(strcmp('FluenceRegion', {sHead.Atlas.Scouts.Label}));
if ~isempty(limiting_scout)
    head_vertices = intersect(head_vertices, limiting_scout.Vertices);
end

if save_in_db && ...,
   ~any(strcmp(['From cortical ' cortex_scout.sScout.Label '(' num2str(extent_m*100) ' cm)']...,
    ,{sHead.Atlas.Scouts.Label}))
    scout_idx = size(sHead.Atlas.Scouts,2) + 1;
    sHead.Atlas.Scouts(scout_idx) = db_template('Scout');
    sHead.Atlas.Scouts(scout_idx).Vertices = head_vertices';
    sHead.Atlas.Scouts(scout_idx).Seed = head_vertices(1);
    sHead.Atlas.Scouts(scout_idx).Color = [0,0,0];
    sHead.Atlas.Scouts(scout_idx).Label = ['From cortical ' cortex_scout.sScout.Label ...
                                           '(' num2str(extent_m*100) ' cm)'];
    bst_save(file_fullpath(sSubject.Surface(sSubject.iScalp).FileName), sHead, 'v7');
    db_save();
end

end