function varargout = process_nst_ComputeFluencesforOptodes( varargin )

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


%addpath(genpath('/NAS/home/zhe_cai/zhengchen/Softwares/iso2mesh'));

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute fluences for optodes by MCXlab';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1200;
sProcess.Description = 'http://mcx.sourceforge.net/cgi-bin/index.cgi';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'data', 'raw'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

sProcess.options.help.Comment = ['<B>This process uses native MEX version of Monte Carlo eXtreme (MCX) to solve the fluences of each optode.</B><BR><BR>' ...
                                     '<B>For more details plese refer to Qianqian Fang and David A. Boas, <BR>' ...
                                     '<B>"Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated by Graphics Processing Units".</B> <BR>', ...
                                     '<B>Opt. Express, vol. 17, issue 22, pp. 20178-20190 (2009)</B><BR>', ...
                                     '<B>For technical details please refer to mcx homepage (http://mcx.sourceforge.net/cgi-bin/index.cgi)</B><BR><BR>'];
sProcess.options.help.Type    = 'label';

ref_length = 10;

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
    sProcess.options.outputdir.Comment = 'Output folder for fluence:';
    sProcess.options.outputdir.Type    = 'filename';
    sProcess.options.outputdir.Value   = SelectOptions;

sProcess.options.mcxlab_gpuid.Comment = pad('cfg.gpuid: ', ref_length);
sProcess.options.mcxlab_gpuid.Type = 'value';
sProcess.options.mcxlab_gpuid.Value={2,'',0};

sProcess.options.mcxlab_autopilot.Comment = pad('cfg.autopilot: ', ref_length);
sProcess.options.mcxlab_autopilot.Type = 'value';
sProcess.options.mcxlab_autopilot.Value={1,'',0};

sProcess.options.mcxlab_respin.Comment = pad('cfg.respin: ', ref_length);
sProcess.options.mcxlab_respin.Type = 'value';
sProcess.options.mcxlab_respin.Value={1,'',0};

sProcess.options.mcxlab_nphoton.Comment = pad('cfg.nphoton: ', ref_length);
sProcess.options.mcxlab_nphoton.Type = 'value';
sProcess.options.mcxlab_nphoton.Value={1,'million photons', 0};
% 
sProcess.options.mcxlab_unitinmm.Comment = pad('cfg.unitinmm: ', ref_length);
sProcess.options.mcxlab_unitinmm.Type = 'value';
sProcess.options.mcxlab_unitinmm.Value={1,'',0};   % defines the length unit for a grid ( voxel) edge length [1.0]
% 
sProcess.options.mcxlab_isreflect.Comment = pad('cfg.isreflect: ', ref_length);
sProcess.options.mcxlab_isreflect.Type = 'value';
sProcess.options.mcxlab_isreflect.Value={1,'',0}; % reflection at exterior boundary
% 
sProcess.options.mcxlab_isrefint.Comment = pad('cfg.isrefint: ', ref_length);
sProcess.options.mcxlab_isrefint.Type = 'value';
sProcess.options.mcxlab_isrefint.Value={1,'',0};   % 1-index mismatch at inner boundaries, [0]-matched index

% time-domain simulation parameters
sProcess.options.mcxlab_tstart.Comment = pad('cfg.tstart: ', ref_length);
sProcess.options.mcxlab_tstart.Type = 'value';
sProcess.options.mcxlab_tstart.Value={0,'',0};
% 
sProcess.options.mcxlab_tend.Comment = pad('cfg.tend:', ref_length);
sProcess.options.mcxlab_tend.Type = 'value';
sProcess.options.mcxlab_tend.Value={5,'ns',1};

sProcess.options.mcxlab_tstep.Comment = pad('cfg.tstep: ', ref_length);
sProcess.options.mcxlab_tstep.Type = 'value';
sProcess.options.mcxlab_tstep.Value={5,'ns',1};
% % nornalisation
sProcess.options.mcxlab_isnormalized.Comment = pad('cfg.isnormalized: ', ref_length);
sProcess.options.mcxlab_isnormalized.Type = 'value';
sProcess.options.mcxlab_isnormalized.Value={1,'',0};
% 
% 
sProcess.options.mcxlab_flag_swepgpu.Comment = 'Swep GPU ';
sProcess.options.mcxlab_flag_swepgpu.Type = 'checkbox';
sProcess.options.mcxlab_flag_swepgpu.Value= 0;

sProcess.options.mcxlab_flag_autoOP.Comment = 'Use default OpticalProperties';
sProcess.options.mcxlab_flag_autoOP.Type = 'checkbox';
sProcess.options.mcxlab_flag_autoOP.Value = 1;

sProcess.options.mcxlab_flag_thresh.Comment = 'Set threshold for fluences (reduce file size)';
sProcess.options.mcxlab_flag_thresh.Type = 'checkbox';
sProcess.options.mcxlab_flag_thresh.Value = 1;
sProcess.options.mcxlab_thresh_value.Comment = 'Threshold for fluences';
sProcess.options.mcxlab_thresh_value.Type = 'value';
sProcess.options.mcxlab_thresh_value.Value = {1,'1e-6(1/mm2/s)',3};

end

function s = pad(s, len)
s = ['<PRE>' s strjoin(repmat({'&nbsp;'}, 1, len-length(s)), '') '</PRE>'];
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

%do_export_fluences = sProcess.options.do_export_fluence_vol.Value;

ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
if ~isfield(ChannelMat.Nirs, 'Wavelengths')
    bst_error(['Head model importation works only for dOD data ' ...
        ' (eg do not use MBLL prior to this process)']);
    return;
end

%addpath(genpath(sProcess.options.mcxpath.Value));


%% Retrieve list of vertex indexes corresponding to current montage
% Load channel file
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);

% Load head mesh
[sSubject, iSubject] = bst_get('Subject', sInputs.SubjectName);
head_mesh_fn = sSubject.Surface(sSubject.iScalp).FileName;
sHead = in_tess_bst(head_mesh_fn);

% Load anat mri
sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
% Load segmenation 
iseg = 0;

for ianat = 1:size(sSubject.Anatomy,2)
    if any(strcmp(sSubject.Anatomy(ianat).Comment, 'segmentation_5tissues'))
        iseg = ianat;
    end
end
if iseg == 0
    bst_error('ERROR: Please import segmentation file as MRI and rename it as "segmentation_5tissues"');
end
seg = in_mri_bst(sSubject.Anatomy(iseg).FileName);

% Retrieve optode coordinates
% Load channel file
[pair_names, tt, tt, pair_sd_idx, src_locs, src_ids, src_chans, det_locs, det_ids, det_chans] = ...
    explode_channels(ChannelMat);
nb_sources = size(src_locs, 1);
nb_dets = size(det_locs, 1);
nb_pairs = length(pair_names);
nb_wavelengths = length(ChannelMat.Nirs.Wavelengths);

% Find closest head vertices (for which we have fluence data)
% Put everything in mri referential
head_vertices_mri = cs_convert(sMri, 'scs', 'mri', sHead.Vertices) * 1000;
src_locs_mri = cs_convert(sMri, 'scs', 'mri', src_locs) * 1000;
det_locs_mri = cs_convert(sMri, 'scs', 'mri', det_locs) * 1000;
head_normals = computeVertNorm(head_vertices_mri,sHead.Faces);
src_hvidx = knnsearch(head_vertices_mri, src_locs_mri);
det_hvidx = knnsearch(head_vertices_mri, det_locs_mri);

%TODO: compute projection errors -> if too large, yell at user

%% Compute fluences from mcxlab
%=========================================================================
% mcxlab setup
%=========================================================================
% GPU thread configuration
flag_swepgpu = sProcess.options.mcxlab_flag_swepgpu.Value;
flag_autoOpticalProperties = sProcess.options.mcxlab_flag_autoOP.Value;
flag_thresh_fluences = sProcess.options.mcxlab_flag_thresh.Value;
thresh_value = sProcess.options.mcxlab_thresh_value.Value{1}.*1e-6;

cfg.gpuid = sProcess.options.mcxlab_gpuid.Value{1};
cfg.autopilot = sProcess.options.mcxlab_autopilot.Value{1};
cfg.respin = sProcess.options.mcxlab_respin.Value{1};
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 
cfg.nphoton=sProcess.options.mcxlab_nphoton.Value{1}.*1e6;
cfg.vol=seg.Cube; 
cfg.unitinmm=sProcess.options.mcxlab_unitinmm.Value{1};   % defines the length unit for a grid ( voxel) edge length [1.0]
cfg.isreflect=sProcess.options.mcxlab_isreflect.Value{1}; % reflection at exterior boundary
cfg.isrefint=sProcess.options.mcxlab_isrefint.Value{1};   % 1-index mismatch at inner boundaries, [0]-matched index
% time-domain simulation parameters
cfg.tstart=sProcess.options.mcxlab_tstart.Value{1}.*1e-9;
cfg.tend=sProcess.options.mcxlab_tend.Value{1}.*1e-9;
cfg.tstep=sProcess.options.mcxlab_tstep.Value{1}.*1e-9;
% nornalisation
cfg.isnormalized=sProcess.options.mcxlab_isnormalized.Value{1}; % [1]-normalize the output flux to unitary source

% project optodes position on cfg.vol
options.proj.stepAlongNorm=1;
options.meshes.skin.vertices=head_vertices_mri;
options.meshes.skin.faces=sHead.Faces;
src_pos=mfip_projectPosInVolume(cfg.vol,head_vertices_mri(src_hvidx,:)+1,head_normals(src_hvidx,:),options,'Display',0,'Text',0);
det_pos=mfip_projectPosInVolume(cfg.vol,head_vertices_mri(det_hvidx,:)+1,head_normals(det_hvidx,:),options,'Display',0,'Text',0);
%TODO: set thresh flag in BST for fluences

for isrc = 1:nb_sources
    cfg.issrcfrom0=0; % [0]- first voxel is [1 1 1] ; 1-first voxel is [0 0 0],
    cfg.srcpos=src_pos(isrc,:);    
    cfg.srcdir=head_normals(src_hvidx(isrc),:);
    for iwl=1:nb_wavelengths
        wl = ChannelMat.Nirs.Wavelengths(iwl);
        if flag_autoOpticalProperties
        cfg.prop=mfip_getOpticalProperties_5layers(1,num2str(wl));
        end
        % running simulation
        fprintf('Running Monte Carlo simulation by MCXlab for Source %g ... \n',isrc);
        if flag_swepgpu
        switch cfg.gpuid 
            case 1
               cfg.gpuid = 2;
            case 2
               cfg.gpuid = 1;
        end
        end
        fluenceRate = mcxlab(cfg); % fluence rate [1/mm2 s]
        fluence.data=fluenceRate.data.*cfg.tstep;  clear fluenceRate
        if flag_thresh_fluences
           fluence.data(fluence.data<thresh_value) = 0;
        end
        volDim=size(fluence.data);      
        fluence_vol = fluence.data;
        fluence_flat_sparse_vol = sparse(double(fluence_vol(:)));
        reference_voxel_index = cfg.srcpos;
        fluence_fn = process_nst_import_head_model('get_fluence_fn', src_hvidx(isrc), wl);
        output_dir = sProcess.options.outputdir.Value{1};
        save([output_dir,'/',fluence_fn],'fluence_flat_sparse_vol','reference_voxel_index');        
    end
end

for idet = 1:nb_dets
    cfg.issrcfrom0=0; % [0]- first voxel is [1 1 1] ; 1-first voxel is [0 0 0],
    cfg.srcpos=det_pos(idet,:);    
    cfg.srcdir=head_normals(det_hvidx(idet),:);    
    for iwl=1:nb_wavelengths
        wl = ChannelMat.Nirs.Wavelengths(iwl);
        if flag_autoOpticalProperties
        cfg.prop=mfip_getOpticalProperties_5layers(1,num2str(wl));
        end
        % running simulation
        fprintf('Running Monte Carlo simulation by MCXlab for Detector %g ... \n',idet);
        if flag_swepgpu
        switch cfg.gpuid 
            case 1
               cfg.gpuid = 2;
            case 2
               cfg.gpuid = 1;
        end
        end
        fluenceRate = mcxlab(cfg); % fluence rate [1/mm2 s]
        fluence.data=fluenceRate.data.*cfg.tstep;  clear fluenceRate
        if flag_thresh_fluences
           fluence.data(fluence.data<thresh_value) = 0;
        end
        volDim=size(fluence.data);
        fluence_vol = fluence.data;
        fluence_flat_sparse_vol = sparse(double(fluence_vol(:)));
        reference_voxel_index = cfg.srcpos;
        fluence_fn = process_nst_import_head_model('get_fluence_fn', det_hvidx(idet), wl);
        output_dir = sProcess.options.outputdir.Value{1};
        save([output_dir,'/',fluence_fn],'fluence_flat_sparse_vol','reference_voxel_index');
    end
end
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


function [pair_names, pair_loc, pair_ichans, pair_sd_indexes, ...
    src_coords, src_ids, src_ichans, ...
    det_coords, det_ids, det_ichans] = explode_channels(channel_def)
%% Explode channel data according to pairs, sources and detectors
% Args
%    - channel_def: struct
%        Definition of channels as given by brainstorm
%        Used fields: Channel
%
% TOCHECK WARNING: uses containers.Map which is available with matlab > v2008
%
%  Outputs:
%     - pair_names: cell array of str, size: nb_pairs
%         Pair names, format: SXDX
%     - pair_loc: array of double, size: nb_pairs x 3 x 2
%         Pair localization (coordinates of source and detector)
%     - pair_ichans: matrix of double, size: nb_pairs x nb_wavelengths
%         Input channel indexes grouped by pairs
%     - pair_sd_indexes: matrix of double, size: nb_pairs x 2
%         1-based continuours indexes of sources and detectors for each
%         sources.
%     - src_coords:   nb_sources x 3
%         Source coordinates, indexed by 1-based continuous index
%         To access via source ID, as read from pair name:
%             src_coords(src_id2idx(src_ID),:)
%     - src_ids: 1d array of double, size: nb_sources
%         vector of source ids (as used in pair name)
%     - src_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the source belongs (indexed by 1-based
%         continuous index).
%     - det_coords:   nb_detectors x 3
%         Detector coordinates, indexed by 1-based continuous index
%         To access via detector ID, as used in pair name:
%             det_coords(det_id2idx(det_ID),:)
%     - det_ids: 1d array of double, size: max_detector_id (hashing vector)
%         vector of detector ids (as used in pair name)
%     - det_chans: cellarray of 1d array of double, size: nb_sources
%         Channel indexes to which the detector belongs (indexed by 1-based
%         continuous index).

MT_OD = 1;
MT_HB = 2;

if isfield(channel_def.Nirs, 'Wavelengths')
    nb_measures = length(channel_def.Nirs.Wavelengths);
    measure_type = MT_OD;
else
    nb_measures = length(channel_def.Nirs.Hb);
    measure_type = MT_HB;
end

pair_to_chans = containers.Map();
pair_to_sd = containers.Map();
src_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
src_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_to_chans = containers.Map('KeyType', 'double', 'ValueType', 'any');
det_coords_map = containers.Map('KeyType', 'double', 'ValueType', 'any');
for ichan=1:length(channel_def.Channel)
    if strcmp(channel_def.Channel(ichan).Type, 'NIRS')
        chan_name = channel_def.Channel(ichan).Name;
        if measure_type == MT_OD
            iwl = strfind(chan_name, 'WL');
            pair_name = chan_name(1:iwl-1);
            wl = str2double(chan_name(iwl+2:end));
            imeasure = channel_def.Nirs.Wavelengths==wl;
        else
            ihb = strfind(chan_name, 'Hb');
            pair_name = chan_name(1:ihb-1);
            imeasure = strcmp(chan_name(ihb:end), channel_def.Nirs.Hb);
        end
        
        if pair_to_chans.isKey(pair_name)
            measures = pair_to_chans(pair_name);
        else
            measures = zeros(1, nb_measures);
        end
        measures(imeasure) = ichan;
        pair_to_chans(pair_name) = measures;
        
        
        [src_id, det_id] = split_pair_name(pair_name);
        pair_to_sd(pair_name) = [src_id, det_id];
        if src_to_chans.isKey(src_id)
            src_to_chans(src_id) = [src_to_chans(src_id) ichan];
        else
            src_to_chans(src_id) = ichan;
            src_coords_map(src_id) = channel_def.Channel(ichan).Loc(:, 1);
        end
        if det_to_chans.isKey(det_id)
            det_to_chans(det_id) = [det_to_chans(det_id) ichan];
        else
            det_to_chans(det_id) = ichan;
            det_coords_map(det_id) = channel_def.Channel(ichan).Loc(:, 2);
        end
        
    end
end

src_coords = cell2mat(src_coords_map.values)';
src_ichans = src_to_chans.values;
src_ids = cell2mat(src_coords_map.keys);

det_coords = cell2mat(det_coords_map.values)';
det_ichans = det_to_chans.values;
det_ids = cell2mat(det_coords_map.keys);

nb_pairs = pair_to_chans.size(1);
pair_names = pair_to_chans.keys;
pair_ichans = zeros(nb_pairs, nb_measures);
pair_loc = zeros(nb_pairs, 3, 2);
pair_sd_indexes = zeros(nb_pairs, 2);
for ipair=1:nb_pairs
    p_indexes = pair_to_chans(pair_names{ipair});
    pair_ichans(ipair, :) = p_indexes;
    pair_loc(ipair, : , :) = channel_def.Channel(pair_ichans(ipair, 1)).Loc;
    sdi = pair_to_sd(pair_names{ipair});
    pair_sd_indexes(ipair, 1) = find(src_ids==sdi(1));
    pair_sd_indexes(ipair, 2) = find(det_ids==sdi(2));
end

end

function [isrc, idet] = split_pair_name(pair_name)
pair_re = 'S([0-9]{1,2})D([0-9]{1,2})';
toks = regexp(pair_name, pair_re , 'tokens');
isrc = str2double(toks{1}{1});
idet = str2double(toks{1}{2});
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


%==========================================================================
% display 
%==========================================================================
if flag_display
    FV=opt.meshes.skin;
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


function fn = protect_fn_str(sfn)
fn = strrep(sfn, ' ', '_');
fn = strrep(fn, '"', '');
fn = strrep(fn, ':', '_');
fn = strrep(fn, '(', '_');
fn = strrep(fn, ')', '_');
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



