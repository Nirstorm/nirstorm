function varargout = process_nst_cpt_fluences( varargin )

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
% Authors: Edouard Delaire (2021)
% Thomas Vincent, ZhengChen Cai (2017)
%
%
eval(macro_method);
end

function sProcess = GetDescription() 
% Description the process
sProcess.Comment     = 'Compute fluences';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1402;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import','data', 'raw'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'data', 'raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;
sProcess.isSeparator = 0;

sProcess.options.subjectname.Comment = 'Subject name:';
sProcess.options.subjectname.Type    = 'subjectname';
sProcess.options.subjectname.Value   = '';


sProcess.options.fluencesCond.Comment = {'panel_nst_fluences', 'Fluence options'};
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
function Comment = FormatComment(sProcess) 
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
%isOk = bst_plugin('Load','mcxlab');
[isOk,] = bst_plugin('Load', sProcess.options.fluencesCond.Value.software);
if ~isOk 
    OutputFiles = {};
    bst_error(['Unable to load ' sProcess.options.fluencesCond.Value.software ]);
    return;
end

% Get scout vertices & load head mesh
if strcmp(sProcess.options.fluencesCond.Value.surface,'montage')
    SubjectName = sProcess.options.fluencesCond.Value.SubjectName;
else
    SubjectName = sProcess.options.subjectname.Value;
end    

sSubject = bst_get('Subject',SubjectName);
if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
    bst_error('No available Cortex and Head surface for this subject.');
    return;
end    

sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);


if strcmp(sProcess.options.fluencesCond.Value.surface,'cortex')
    i_atlas = strcmp({sCortex.Atlas.Name},sProcess.options.fluencesCond.Value.Atlas);
    cortex_to_scalp_extent = sProcess.options.fluencesCond.Value.Extent;

    rois = strsplit(sProcess.options.fluencesCond.Value.ROI,',');
    head_vertices = [];
    
    for i_roi = 1:length(rois)
        i_scout = strcmp({sCortex.Atlas(i_atlas).Scouts.Label}, strtrim(rois(i_roi)));

        cortex_scout            = struct();
        cortex_scout.sSubject   = sSubject;
        cortex_scout.sScout     = sCortex.Atlas(i_atlas).Scouts(i_scout);

        head_vertices = union(head_vertices,  ...
                        bst_getoutvar(1,@proj_cortex_scout_to_scalp,cortex_scout,cortex_to_scalp_extent.*0.01, 1));
    end
elseif strcmp(sProcess.options.fluencesCond.Value.surface,'head')
    i_atlas = strcmp({sHead.Atlas.Name},sProcess.options.fluencesCond.Value.Atlas);
      
    rois = strsplit(sProcess.options.fluencesCond.Value.ROI,',');
    head_vertices = [];
    
    for i_roi = 1:length(rois)
        i_scout = strcmp({sHead.Atlas(i_atlas).Scouts.Label},  strtrim(rois(i_roi)));
        head_vertices = union(head_vertices, sHead.Atlas(i_atlas).Scouts(i_scout).Vertices);    
    end  
    
    iHeadAtlas      = find(strcmp({sHead.Atlas.Name}, 'User scouts'));
    iScoutExclude   = find(strcmp({sHead.Atlas(iHeadAtlas).Scouts.Label},'FluenceExclude'));
    exclude_scout   = sHead.Atlas(iHeadAtlas).Scouts(iScoutExclude);

    if ~isempty(exclude_scout)
        head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
    end

else
    
    ChannelFile = sProcess.options.fluencesCond.Value.ChannelFile;
    ChannelMat = in_bst_channel(ChannelFile);
    
    if ~isfield(ChannelMat.Nirs, 'Wavelengths')
        bst_error(['Fluence can only be computed on raw or dOD data ' ...
            ' (eg do not use MBLL prior to this process)']);
        return;
    end
    sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
    [head_vertices] = proj_montage_to_scalp(ChannelMat, sHead, sMri);
    
end

    
OutputFiles = Compute(sProcess, sSubject, sHead, head_vertices);
end

function [head_vertices] = proj_montage_to_scalp(ChannelMat, sHead, sMri)
    % Retrieve optode coordinates
    % Load channel file
    montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel);
    src_locs = montage_info.src_pos;
    det_locs = montage_info.det_pos;

    [src_hvidx, det_hvidx] = process_nst_import_head_model('get_head_vertices_closest_to_optodes',...
                                                    sMri, sHead, src_locs, det_locs);

    head_vertices = [ src_hvidx' , det_hvidx'];
end

function [head_vertices, sHead, sSubject] = proj_cortex_scout_to_scalp(cortex_scout, extent_m, save_in_db)

if nargin < 3
    save_in_db = 0;
end
sSubject = cortex_scout.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
dis2head = nst_pdist(sHead.Vertices, sCortex.Vertices(cortex_scout.sScout.Vertices,:));
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

function OutputFiles = Compute(sProcess, sSubject, sHead, valid_vertices)

OutputFiles = {};
options = sProcess.options.fluencesCond.Value;

% Load anat mri
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

% Load segmentation 
segmentation_name = 'segmentation_5tissues';
iseg = strcmp({sSubject.Anatomy.Comment}, segmentation_name);
if ~any(iseg)
    bst_error(sprintf('ERROR: given segmentation "%s" not found for subject "%s"', ...
                      segmentation_name, sSubject.Name));
    return
end

sSegmentation   = in_mri_bst(sSubject.Anatomy(iseg).FileName);
if ~isequal(sMri.Voxsize,sSegmentation.Voxsize)
    bst_report('Error', sProcess, [], 'MRI and Segmentation have different voxel size');
    return
end

if ~isfield(sSegmentation,'Labels') || isempty(sSegmentation.Labels)
    bst_report('Error', sProcess, [], ['BST> Invalid atlas "' segmentation_name '": does not contain any labels.']);
    return;
end
tissues         = sSegmentation.Labels;

% Find closest head vertices (for which we have fluence data)
% Put everything in mri referential
head_vertices_mri   = cs_convert(sMri, 'scs', 'voxel', sHead.Vertices);
head_normals        = tess_normals(head_vertices_mri,sHead.Faces); %Use brainstorm
head_normals        = -head_normals;

%load wavelength info 
wavelengths_input = options.wavelengths;
try
    scan_res = textscan(wavelengths_input, '%d,');
catch
    bst_report('Error', sProcess, [], 'List of wavelengths must be integers separated by comas');
    return
end

wavelengths     = double(scan_res{1}');
nb_wavelengths  = length(wavelengths);
nb_vertex       = length(valid_vertices);

%% Compute fluences from mcxlab
%=========================================================================
% mcxlab setup
%=========================================================================

flag_overwrite_fluences = options.mcxlab_overwrite_fluences;

cfg.gpuid       = options.mcxlab_gpuid;
cfg.autopilot   = 1;
cfg.respin      = 1;
% set seed to make the simulation repeatible
cfg.seed        = hex2dec('623F9A9E'); 
cfg.nphoton     = options.mcxlab_nphoton*1e6;
cfg.vol         = sSegmentation.Cube; % segmentation
cfg.unitinmm    = sSegmentation.Voxsize(1);   % defines the length unit for a grid ( voxel) edge length [1.0]
cfg.isreflect   = 1; % reflection at exterior boundary
cfg.isrefint    = 1; % 1-index mismatch at inner boundaries, [0]-matched index
% time-domain simulation parameters
cfg.tstart      = 0;
cfg.tend        = 5e-9;
cfg.tstep       = 5e-9;
% nornalisation
cfg.isnormalized= 1; % [1]-normalize the output flux to unitary source

cfg.issrcfrom0  = 0; % [0]- first voxel is [1 1 1] ; 1-first voxel is [0 0 0],
% project optodes position on cfg.vol
options.proj.stepAlongNorm      = 1;
options.meshes.skin.vertices    = head_vertices_mri;
options.meshes.skin.faces       = sHead.Faces;

[vertex_pos,invalid_id]=mfip_projectPosInVolume(cfg.vol,head_vertices_mri(valid_vertices,:)+1,head_normals(valid_vertices,:),options,'Display',0,'Text',0);

if ~isempty(invalid_id)
    
    bst_report('Warning', sProcess, [], 'Some vertices could not be project into the volume and were discarded (check scout Wrong vertex)');

    scout_idx = size(sHead.Atlas.Scouts,2) + 1;
    sHead.Atlas.Scouts(scout_idx) = db_template('Scout');
    sHead.Atlas.Scouts(scout_idx).Vertices = valid_vertices(invalid_id)';
    sHead.Atlas.Scouts(scout_idx).Seed = sHead.Atlas.Scouts(scout_idx).Vertices(1);
    sHead.Atlas.Scouts(scout_idx).Color = [0,0,0];
    sHead.Atlas.Scouts(scout_idx).Label = sprintf('Wrong vertex for ROI %s', options.ROI);

    bst_save(file_fullpath(sSubject.Surface(sSubject.iScalp).FileName), sHead, 'v7');
    db_save();
    % remove vertices from fluence

    valid_vertices = setdiff(valid_vertices,valid_vertices(invalid_id));
end
    
tic
bst_progress('start', 'Compute fluences', sprintf('Computing fluences for %d vertices and %d wavelengths', nb_vertex, nb_wavelengths), ...
             1, nb_vertex * nb_wavelengths);
for ivertx = 1:nb_vertex
    cfg.srcpos=vertex_pos(ivertx,:);    
    cfg.srcdir=head_normals(valid_vertices(ivertx),:);
    for iwl=1:nb_wavelengths
        wl = wavelengths(iwl);
        fluence_fn = process_nst_import_head_model('get_fluence_fn', valid_vertices(ivertx), wl);
        output_dir = options.outputdir;
        if exist(fullfile(output_dir, fluence_fn), 'file') == 2 && flag_overwrite_fluences == 0
            fprintf('Fluence for head vertex %g exists in the target folder\n',valid_vertices(ivertx));
            bst_progress('inc', 1);
        else

            cfg.prop = nst_get_tissues_optical_properties(tissues,wl);
            fprintf('Running Monte Carlo simulation by MCXlab for head vertex %g ... \n',valid_vertices(ivertx));
            
            if strcmp(options.software, 'mcxlab-cuda')
                fluenceRate = mcxlab(cfg); % fluence rate [1/mm2 s]
            else
                fluenceRate = mcxlabcl(cfg); % fluence rate [1/mm2 s]
            end
            fluenceRate.data=fluenceRate.data.*cfg.tstep; 
            
            if isempty(fluenceRate.data) 
                bst_error(sprintf('ERROR:Fluence %d cannot be computed (see command windows)',valid_vertices(ivertx)));
                return;
            elseif isinf(fluenceRate.stat.normalizer) && strcmp(bst_get('OsType'),'mac64arm')
                bst_error('ERROR: MCXlabCL reach time limit (Try again with less photons)');
                return;
            elseif isinf(fluenceRate.stat.normalizer)
                bst_error(sprintf('ERROR:Fluence %d cannot be computed (see command windows)',valid_vertices(ivertx)));
                return;
            end
            
            fluence_vol = fluenceRate.data;

            fluence_flat_sparse_vol = sparse(double(fluence_vol(:)));
            reference_voxel_index   = cfg.srcpos;

            save(fullfile(output_dir, fluence_fn), 'fluence_flat_sparse_vol','reference_voxel_index');
            bst_progress('inc', 1);
        end
    end
end
bst_progress('stop', 1);
disp('');
disp([num2str(nb_vertex * nb_wavelengths) ' fluence volumes computed in ' num2str(toc) ' seconds']);

end

function [normalV] = computeVertNorm(Vertices,Faces)
normalF = surfacenorm(Vertices,Faces);
nvert = size(Vertices,1);
nface = size(Faces,1);
normalV = zeros(nvert,3);
for i=1:nface
    F = Faces(i,:);
    for j=1:3
        normalV(F(j),:) = normalV(F(j),:) + normalF(i,:);
    end
end
% normalize
d = sqrt( sum(normalV.^2,2) ); d(d<eps)=1;
normalV = normalV ./ repmat( d, 1,3 );
% put normals inward
v = Vertices - repmat(mean(Vertices,2), 1,3);
s = sum( v.*normalV, 2 );
if sum(s>0)>sum(s<0)
    % flip
    normalV = -normalV;
end
end

function [pos, invalid_id] = mfip_projectPosInVolume(vol,pos,normals,options,varargin)

%==========================================================================
% parse optional parameters
%==========================================================================
p = inputParser; 

% add Optional
p.addOptional ('options',struct,@(x)1)

% param/value
p.addParamValue ('Display',0,@(x)1)
p.addParamValue ('Text',0,@(x)1)

% parse
p.parse(options,varargin{:})

opt=p.Results.options;

flag_display=p.Results.Display;
flag_text=p.Results.Text;


% parse optional parameters
% skin
if isempty(opt.meshes.skin)
flag_display=0;
end
% stepAlongNorm
stepAlongNorm=opt.proj.stepAlongNorm;
FV=opt.meshes.skin;

invalid_id = [];
%==========================================================================
% display 
%==========================================================================
if flag_display
    %FV=opt.meshes.skin
    

    hFig1=figure('name','before projection');
    hPatches = patch(FV,'Parent',gca);
    set(hPatches,...
        'facecolor', [.5 .5 .5], ...
        'EdgeColor', 'none',...
        'FaceAlpha', 0.3);
    axis equal;
    mfip_setAxesOrientation(gca,'RAS')
    hold on
    plot3(pos(:,1),pos(:,2),pos(:,3),'.r','MarkerSize',25)
    quiver3(pos(:,1),pos(:,2),pos(:,3),...
    10*normals(:,1),10*normals(:,2),10*normals(:,3))
    if  flag_text
        for ii=1:size(pos,1)
            symbol = sprintf('P%g',ii);
            text(pos(ii,1),pos(ii,2),pos(ii,3),symbol,'FontSize',14);
        end
    end
end

%==========================================================================
% project optodes on non background tissue :
% if the specified optode position is outside the domain, move the optode
% along the initial vector until it hits the domain
%==========================================================================
[dim1, dim2, dim3]=size(vol);


pos=round(pos); % if the voxel coordinates are not integers change them for the nearest integer
pos_orig=pos;
idx = pos<1; pos(idx)=1; clear idx
idx = pos(:,1)>dim1; pos(idx,1)=dim1; clear idx
idx = pos(:,2)>dim2; pos(idx,2)=dim2; clear idx
idx = pos(:,3)>dim3; pos(idx,3)=dim3; clear idx

%==========================================================================
% searching non-zero voxel failed along the incident vector
%==========================================================================

% find linear index in the volume matrix
hFig1 = [];

for iOpt=1:size(pos,1)
    if vol(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3))==0 % 0 is for background in fuzzy vol
        fprintf('opt %g (%f %f %f) is located outside the domain\n',iOpt,pos(iOpt,:));
       ite=0;
        while vol(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3))==0 & ite<400;
            ite=ite+1;
            pos(iOpt,:)=pos(iOpt,:)+stepAlongNorm*normals(iOpt,:);
            pos=round(pos); 
            idx = pos<1; pos(idx)=1; clear idx
            idx = pos(:,1)>dim1; pos(idx,1)=dim1; clear idx
            idx = pos(:,2)>dim2; pos(idx,2)=dim2; clear idx
            idx = pos(:,3)>dim3; pos(idx,3)=dim3; clear idx
            vol(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3));
        end
        
        % =================================================================
        % Controlling bad optodes
        % =================================================================
        
        if  vol(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3))==0
            fprintf('opt %g (%f %f %f) cannot be fixed ; stopping\n',iOpt,pos(iOpt,:));
%             [hFig1, iFig, iDS] = bst_figures('GetFigure', hFig1);
%             if isempty(hFig1)
%                 hFig1=figure('name','bad optode');
%                 
%             subplot(1,1,1)
%             hPatches = patch(FV,'Parent',gca);
%             set(hPatches,...
%                 'facecolor', [.5 .5 .5], ...
%                 'EdgeColor', 'none',...
%                 'FaceAlpha', 0.3);
%             axis equal;
%             mfip_setAxesOrientation(gca,'RAS')
%             end
%             hold on
%             plot3(pos_orig(iOpt,1),pos(iOpt,2),pos(iOpt,3),'.b','MarkerSize',25)
%             plot3(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3),'.r','MarkerSize',25)
%             quiver3(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3),10*normals(iOpt,1),10*normals(iOpt,2),10*normals(iOpt,3))
%             
%           
%              x=1:dim1; y=1:dim2; z=1:dim3;
%              [a b c]=ndgrid(1:dim1,1:dim2,1:dim3);
%              volVoxCoord=[a(:) b(:) c(:)]; clear a b c
%              
%              viewMode='sagital';
%              switch viewMode
%                  case 'axial'
%              [mx,my] = meshgrid(x,y);
%              coeff(1)=0;
%              coeff(2)=0;
%              coeff(3)=1;
%              coeff(4)=-pos_orig(iOpt,3);
%              imData=squeeze(vol(:,:,pos_orig(iOpt,3)))';
%              mz = (coeff(1)*mx+ coeff(2)*my + coeff(4)) / -coeff(3);
%              mz=round(mz);  % the slice will be defined with the nearest voxel of mz
%               case 'sagital'
%              [my,mz] = meshgrid(y,z);
%              coeff(1)=1;
%              coeff(2)=0;
%              coeff(3)=0;
%              coeff(4)=-pos_orig(iOpt,1);
%              imData=squeeze(vol(pos_orig(iOpt,1),:,:))';
%              mx = (coeff(3)*mz+ coeff(2)*my + coeff(4)) / -coeff(1);
%              mx=round(mx);  % the slice will be defined with the nearest voxel of mz
%              case 'coronal'
%              [mx,mz] = meshgrid(x,z);
%              coeff(1)=0;
%              coeff(2)=1;
%              coeff(3)=0;
%              coeff(4)=-pos_orig(iOpt,2);
%              imData=squeeze(vol(:,:,pos_orig(iOpt,2)))';  
%              my = (coeff(3)*mz+ coeff(1)*mx + coeff(4)) / -coeff(2);
%              my=round(my);  % the slice will be defined with the nearest voxel of mz
%              end
%              hold on
%              hsp =surf(mx,my,mz,double(imData),'edgecolor','none');     
%              clear x y z
             if nargout < 2 
                error('Stopped because the algorithm cannot project one optode')
             else
                 invalid_id(end+1) = iOpt;
             end    
        else
        fprintf('fixing position to (%f %f %f)\n\n',pos(iOpt,:));
        end
    end
end

%==========================================================================
% display 
%==========================================================================
if flag_display
    hFig1=figure('name','after projection');
    hPatches = patch(FV,'Parent',gca);
    set(hPatches,...
        'facecolor', [.5 .5 .5], ...
        'EdgeColor', 'none',...
        'FaceAlpha', 0.3);
    axis equal;
    mfip_setAxesOrientation(gca,'RAS')
    hold on
    plot3(pos(:,1),pos(:,2),pos(:,3),'.r','MarkerSize',25)
    
    if  flag_text
        for ii=1:size(pos,1)
            symbol = sprintf('P%g',ii);
            text(pos(ii,1),pos(ii,2),pos(ii,3),symbol,'FontSize',14);
        end
    end
end
end



function mfip_setAxesOrientation (hAxes, orientation)

for iA=1:length(hAxes)
    switch orientation 
            case 'RAS'
        set(get(hAxes(iA),'XLabel'),'String','x:L-->R')
        set(get(hAxes(iA),'YLabel'),'String','y:P-->A')
        set(get(hAxes(iA),'ZLabel'),'String','z:I-->S') 
        set(hAxes(iA),'Xdir','reverse')
        set(hAxes(iA),'Ydir','reverse')
        set(hAxes(iA),'Zdir','normal')

            case 'RPI'
        set(get(hAxes(iA),'XLabel'),'String','x:L-->R')
        set(get(hAxes(iA),'YLabel'),'String','y:A-->P')
        set(get(hAxes(iA),'ZLabel'),'String','z:S-->I') 
        set(hAxes(iA),'Xdir','reverse')
        set(hAxes(iA),'Ydir','normal')
        set(hAxes(iA),'Zdir','reverse')
    end 

end

% non stretch to fill behavior when rotating http://www.mathworks.com/matlabcentral/newsreader/view_thread/9714
axis(hAxes,'vis3d')
axis(hAxes, 'equal')
end