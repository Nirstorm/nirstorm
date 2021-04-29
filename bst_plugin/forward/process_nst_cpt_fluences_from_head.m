function varargout = process_nst_cpt_fluences_from_head( varargin )

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
%
% Main process to compute fluences. All other processes call this one
% via the Compute function: 
%     - process_nst_compute_fluences_from_cortex_scout.m
%     - process_nst_compute_fluences_from_montage.m
%
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute fluences from head scout';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Sources'};
sProcess.Index       = 1404;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

sProcess.options = nst_add_scout_sel_options(struct(), 'head', str_pad('Head scout (search space):',40), ...
                                             'scalp', {'User scouts'}, 1); 

sProcess.options = append_mcxlab_options(sProcess.options);
end

function options = append_mcxlab_options(options)

options.segmentation.Comment = str_pad('Segmentation (Volume)',40);
options.segmentation.Type    = 'text';
options.segmentation.Value = 'segmentation';

options.segmentation_label.Type    = 'radio_line';
options.segmentation_label.Comment   = {'1:skin, 2:skull, 3:CSF, 4:GM, 5:WM', '5: skin,  4: skull, 3: CSF, 2: GM, 1: WM','Segmentation label: '};
options.segmentation_label.Value   = 1;   

options.wavelengths.Comment = str_pad('Wavelengths (nm) [coma-separated list]',40);
options.wavelengths.Type    = 'text';
options.wavelengths.Value = '';

options.help.Comment = ['<B>This process uses native MEX version of Monte Carlo eXtreme (MCX) to solve the fluences of each optode.</B><BR><BR>' ...
    '<B>For more details plese refer to Qianqian Fang and David A. Boas, <BR>' ...
    '<B>"Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated by Graphics Processing Units".</B> <BR>', ...
    '<B>Opt. Express, vol. 17, issue 22, pp. 20178-20190 (2009)</B><BR>', ...
    '<B>For technical details please refer to mcx homepage (http://mcx.sourceforge.net/cgi-bin/index.cgi)</B><BR><BR>'];
options.help.Type    = 'label';

% ref_length = 26;

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
options.outputdir.Comment = 'Output folder for fluence:';
options.outputdir.Type    = 'filename';
options.outputdir.Value   = SelectOptions;

options.mcxlab_gpuid.Comment = str_pad('cfg.gpuid: ',20);
options.mcxlab_gpuid.Type = 'value';
options.mcxlab_gpuid.Value={1,'',0};

options.mcxlab_nphoton.Comment = str_pad('cfg.nphoton: ',20);
options.mcxlab_nphoton.Type = 'value';
options.mcxlab_nphoton.Value={1,'million photons', 0};

options.mcxlab_flag_autoOP.Comment = str_pad('Use default OpticalProperties',20);
options.mcxlab_flag_autoOP.Type = 'checkbox';
options.mcxlab_flag_autoOP.Value = 1;

options.mcxlab_flag_thresh.Comment = str_pad('Set threshold for fluences (reduce file size)',20);
options.mcxlab_flag_thresh.Type = 'checkbox';
options.mcxlab_flag_thresh.Value = 0;

options.mcxlab_thresh_value.Comment = str_pad('Threshold for fluences',20);
options.mcxlab_thresh_value.Type = 'value';
options.mcxlab_thresh_value.Value = {1,'1e-6(1/mm2/s)',0};

options.mcxlab_overwrite_fluences.Comment = str_pad('Overwirte existing fluences',20);
options.mcxlab_overwrite_fluences.Type = 'checkbox';
options.mcxlab_overwrite_fluences.Value = 0;
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
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% Get scout vertices & load head mesh
head_scout_selection = nst_get_option_selected_scout(sProcess.options, 'head');
head_vertices = head_scout_selection.sScout.Vertices;

sSubject = head_scout_selection.sSubject;
sHead = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);

OutputFiles = Compute(sProcess, sSubject, sHead, head_vertices);
end

function OutputFiles = Compute(sProcess, sSubject, sHead, valid_vertices)

OutputFiles = {};

% Load anat mri
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);

% Load segmentation 
iseg = strcmp(sProcess.options.segmentation.Value, {sSubject.Anatomy.Comment});
if ~any(iseg)
    bst_error(sprintf('ERROR: given segmentation "%s" not found for subject "%s"', ...
                      sProcess.options.segmentation.Value, sSubject.Name));
    return
end
seg = in_mri_bst(sSubject.Anatomy(iseg).FileName);
if sProcess.options.segmentation_label.Value == 1
    seg.Cube = nst_prepare_segmentation(seg.Cube,{1,2,3,4,5});
elseif sProcess.options.segmentation_label.Value == 2
    seg.Cube = nst_prepare_segmentation(seg.Cube,{5,4,3,2,1});
end    
% Find closest head vertices (for which we have fluence data)
% Put everything in mri referential
head_vertices_mri = cs_convert(sMri, 'scs', 'mri', sHead.Vertices) * 1000;
% head_normals = computeVertNorm(head_vertices_mri,sHead.Faces); %Use iso2mesh
head_normals = tess_normals(head_vertices_mri,sHead.Faces); %Use brainstorm
% put normals inward
head_normals = -head_normals;

nb_vertex = length(valid_vertices);

%load wavelength info 
wavelengths_input = sProcess.options.wavelengths.Value;
try
    scan_res = textscan(wavelengths_input, '%d,');
catch
    bst_report('Error', sProcess, [], 'List of wavelengths must be integers separated by comas');
    return
end
wavelengths = double(scan_res{1}');
nb_wavelengths = length(wavelengths);

%% Compute fluences from mcxlab
%=========================================================================
% mcxlab setup
%=========================================================================
% GPU thread sadfsconfiguration

flag_autoOpticalProperties = sProcess.options.mcxlab_flag_autoOP.Value;
flag_thresh_fluences = sProcess.options.mcxlab_flag_thresh.Value;
thresh_value = sProcess.options.mcxlab_thresh_value.Value{1}.*1e-6;
flag_overwrite_fluences = sProcess.options.mcxlab_overwrite_fluences.Value;

cfg.gpuid = sProcess.options.mcxlab_gpuid.Value{1};
cfg.autopilot = 1;
cfg.respin = 1;
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 
cfg.nphoton=sProcess.options.mcxlab_nphoton.Value{1}.*1e6;
cfg.vol=seg.Cube; % segmentation
cfg.unitinmm=1;   % defines the length unit for a grid ( voxel) edge length [1.0]
cfg.isreflect=1; % reflection at exterior boundary
cfg.isrefint=1;   % 1-index mismatch at inner boundaries, [0]-matched index
% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
% nornalisation
cfg.isnormalized=1; % [1]-normalize the output flux to unitary source

cfg.issrcfrom0=0; % [0]- first voxel is [1 1 1] ; 1-first voxel is [0 0 0],
% project optodes position on cfg.vol
options.proj.stepAlongNorm=1;
options.meshes.skin.vertices=head_vertices_mri;
options.meshes.skin.faces=sHead.Faces;
vertex_pos=mfip_projectPosInVolume(cfg.vol,head_vertices_mri(valid_vertices,:)+1,head_normals(valid_vertices,:),options,'Display',0,'Text',0);
%det_pos=mfip_projectPosInVolume(cfg.vol,head_vertices_mri(det_hvidx,:)+1,head_normals(det_hvidx,:),options,'Display',0,'Text',0);
%TODO: set thresh flag in BST for fluences

tic
bst_progress('start', 'Compute fluences', sprintf('Computing fluences for %d vertices and %d wavelengths', nb_vertex, nb_wavelengths), ...
             1, nb_vertex * nb_wavelengths);
for ivertx = 1:nb_vertex
    cfg.srcpos=vertex_pos(ivertx,:);    
    cfg.srcdir=head_normals(valid_vertices(ivertx),:);
    for iwl=1:nb_wavelengths
        wl = wavelengths(iwl);
        fluence_fn = process_nst_import_head_model('get_fluence_fn', valid_vertices(ivertx), wl);
        output_dir = sProcess.options.outputdir.Value{1};
        if exist(fullfile(output_dir, fluence_fn), 'file') == 2 && flag_overwrite_fluences == 0
            fprintf('Fluence for head vertex %g exists in the target folder\n',valid_vertices(ivertx));
            bst_progress('inc', 1);
        else
            if flag_autoOpticalProperties
                cfg.prop=mfip_getOpticalProperties_5layers(1,num2str(wl));
            end
            % running simulation
            fprintf('Running Monte Carlo simulation by MCXlab for head vertex %g ... \n',valid_vertices(ivertx));
            
            fluenceRate = mcxlab(cfg); % fluence rate [1/mm2 s]
            fluence.data=fluenceRate.data.*cfg.tstep;  clear fluenceRate
            if flag_thresh_fluences
                fluence.data(fluence.data<thresh_value) = 0;
            end
            volDim=size(fluence.data);
            fluence_vol = fluence.data;
            fluence_flat_sparse_vol = sparse(double(fluence_vol(:)));
            reference_voxel_index = cfg.srcpos;
            %         fluence_fn = process_nst_import_head_model('get_fluence_fn', valid_vertices(ivertx), wl);
            %         output_dir = sProcess.options.outputdir.Value{1};
            save(fullfile(output_dir, fluence_fn), 'fluence_flat_sparse_vol','reference_voxel_index');
            bst_progress('inc', 1);
        end
    end
end
bst_progress('stop', 1);
disp('');
disp([num2str(nb_vertex * nb_wavelengths) ' fluence volumes computed in ' num2str(toc) ' seconds']);
%save([output_dir,'/','montage_region.mat'],'newScout');



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

function [pos] = mfip_projectPosInVolume(vol,pos,normals,options,varargin)

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

%==========================================================================
% display 
%==========================================================================
if flag_display
    %FV=opt.meshes.skin;
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
[dim1 dim2 dim3]=size(vol);


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
            fprintf('opt %g (%f %f %f) cannot be fixed ; stopping\n');
            hFig1=figure('name','bad optode');
            subplot(1,1,1)
            hPatches = patch(FV,'Parent',gca);
            set(hPatches,...
                'facecolor', [.5 .5 .5], ...
                'EdgeColor', 'none',...
                'FaceAlpha', 0.3);
            axis equal;
            mfip_setAxesOrientation(gca,'RAS')
            hold on
            plot3(pos_orig(iOpt,1),pos(iOpt,2),pos(iOpt,3),'.b','MarkerSize',25)
            plot3(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3),'.r','MarkerSize',25)
            quiver3(pos(iOpt,1),pos(iOpt,2),pos(iOpt,3),10*normals(iOpt,1),10*normals(iOpt,2),10*normals(iOpt,3))
            
          
             x=1:dim1; y=1:dim2; z=1:dim3;
             [a b c]=ndgrid(1:dim1,1:dim2,1:dim3);
             volVoxCoord=[a(:) b(:) c(:)]; clear a b c
             
             viewMode='sagital';
             switch viewMode
                 case 'axial'
             [mx,my] = meshgrid(x,y);
             coeff(1)=0;
             coeff(2)=0;
             coeff(3)=1;
             coeff(4)=-pos_orig(iOpt,3);
             imData=squeeze(vol(:,:,pos_orig(iOpt,3)))';
             mz = (coeff(1)*mx+ coeff(2)*my + coeff(4)) / -coeff(3);
             mz=round(mz);  % the slice will be defined with the nearest voxel of mz
              case 'sagital'
             [my,mz] = meshgrid(y,z);
             coeff(1)=1;
             coeff(2)=0;
             coeff(3)=0;
             coeff(4)=-pos_orig(iOpt,1);
             imData=squeeze(vol(pos_orig(iOpt,1),:,:))';
             mx = (coeff(3)*mz+ coeff(2)*my + coeff(4)) / -coeff(1);
             mx=round(mx);  % the slice will be defined with the nearest voxel of mz
             case 'coronal'
             [mx,mz] = meshgrid(x,z);
             coeff(1)=0;
             coeff(2)=1;
             coeff(3)=0;
             coeff(4)=-pos_orig(iOpt,2);
             imData=squeeze(vol(:,:,pos_orig(iOpt,2)))';  
             my = (coeff(3)*mz+ coeff(1)*mx + coeff(4)) / -coeff(2);
             my=round(my);  % the slice will be defined with the nearest voxel of mz
             end
             hold on
             hsp =surf(mx,my,mz,double(imData),'edgecolor','none');     
             clear x y z
            error('Stopped because the algorithm cannot project one optode')
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

function [prop] = mfip_getOpticalProperties_5layers(method,wavelength)

%==========================================================================
%                       INITIALISATION: SET DEFAULTS
%==========================================================================
p = inputParser;
addRequired(p,'method',@isscalar);
addRequired(p,'wavelength',@ischar);
parse(p,method,wavelength);



%==========================================================================
%                       GET OPTICAL PROPERTIES
%==========================================================================
switch method
    case 1
        % Strangman_2003: mus was calculated from mus', conversion in mm-1
        % NB: another useful references : Yamada 2009 JBO 14(6) has a set of absorption and scattering coeff at multiple wavelength for the five layers

        
        if strcmp(wavelength,'685')
            % 690 nm : [mua, mus, g, n]
            prop=[0       0       1    1      % medium 0: the environment
                0.0159  10      0.92 1.37   % medium 1: skin
                0.0101  12.5    0.92 1.37   % medium 2: skull
                0.0004  0.125   0.92 1.37   % medium 3: CSF
                0.0178  15.625  0.92 1.37   % medium 4: gray matter
                0.0178  15.625  0.92 1.37]; % medium 5: white matter
            
       elseif strcmp(wavelength,'690')
            % 690 nm : [mua, mus, g, n]
            prop=[0       0       1    1      % medium 0: the environment
                0.0159  10      0.92 1.37   % medium 1: skin
                0.0101  12.5    0.92 1.37   % medium 2: skull
                0.0004  0.125   0.92 1.37   % medium 3: CSF
                0.0178  15.625  0.92 1.37   % medium 4: gray matter
                0.0178  15.625  0.92 1.37]; % medium 5: white matter
            
        elseif strcmp(wavelength,'830')
            % 830 nm : [mua, mus, g, n]
            prop=[0      0      1     1       % medium 0: the environment
                0.0191 8.25   0.92  1.37    % medium 1: skin
                0.0136 10.75  0.92  1.37    % medium 2: skull
                0.0026 0.125   0.92 1.37    % medium 3: CSF
                0.0186 13.875 0.92  1.37    % medium 4: gray matter
                0.0186 13.875 0.92  1.37]; % medium 5: white matter
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



