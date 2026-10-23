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
% Authors: Edouard Delaire (2021-2026)
% Thomas Vincent, ZhengChen Cai (2017)
%
%
    eval(macro_method);
end

function sProcess = GetDescription() 
    % Description the process
    sProcess.Comment     = 'Compute fluences';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Forward model'};
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
        sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
        
        ChannelFile = unique({sInputs.ChannelFile});
        head_vertices = [];
    
        for iChannel = 1:length(ChannelFile)
    
            ChannelMat = in_bst_channel(ChannelFile{iChannel});
    
            if ~isfield(ChannelMat.Nirs, 'Wavelengths')
                bst_error(['Fluence can only be computed on raw or dOD data ' ...
                    ' (eg do not use MBLL prior to this process)']);
                return;
            end
            
            % Ensure that all the nirs sensors are on the head
            iNirs = good_channel(ChannelMat.Channel, [], 'nirs');
            allLoc = [ChannelMat.Channel(iNirs).Loc]';
            ChanLocPrj = channel_project_scalp(sHead.Vertices, allLoc);
            if any(sqrt(sum( (ChanLocPrj - allLoc).^2, 2)) > 0.001)
                bst_error(['One or more nirs sensors are not placed on the scalp.' 10 ...
                    'Project the sensors to the scalp surface before computing the head model.']);
                return
            end
    
            head_vertices = union(head_vertices,  ...
                                  bst_getoutvar(1, @proj_montage_to_scalp, ChannelMat, sHead, sMri));
    
        end
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
    
    iHeadAtlas    = find(strcmp({sHead.Atlas.Name}, 'User scouts'));
    exclude_scout = sHead.Atlas(iHeadAtlas).Scouts(strcmp('FluenceExclude', {sHead.Atlas(iHeadAtlas).Scouts.Label}));
    if ~isempty(exclude_scout)
        head_vertices = setdiff(head_vertices, exclude_scout.Vertices);
    end
    
    limiting_scout = sHead.Atlas(iHeadAtlas).Scouts(strcmp('FluenceRegion', {sHead.Atlas(iHeadAtlas).Scouts.Label}));
    if ~isempty(limiting_scout)
        head_vertices = intersect(head_vertices, limiting_scout.Vertices);
    end
    
    if save_in_db && ...,
       ~any(strcmp(['From cortical ' cortex_scout.sScout.Label '(' num2str(extent_m*100) ' cm)']...,
        ,{sHead.Atlas(iHeadAtlas).Scouts.Label}))
        scout_idx = size(sHead.Atlas(iHeadAtlas).Scouts,2) + 1;
        sHead.Atlas(iHeadAtlas).Scouts(scout_idx) = db_template('Scout');
        sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Vertices = head_vertices';
        sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Seed = head_vertices(1);
        sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Color = [0,0,0];
        sHead.Atlas(iHeadAtlas).Scouts(scout_idx).Label = ['From cortical ' cortex_scout.sScout.Label ...
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
    if ~all(round(sMri.Voxsize(1:3) .* 1000) == round(sSegmentation.Voxsize(1:3) .* 1000))
        bst_report('Error', sProcess, [], 'MRI and Segmentation have different voxel size');
        return
    end
    
    if ~isfield(sSegmentation,'Labels') || isempty(sSegmentation.Labels)
        bst_report('Error', sProcess, [], ['BST> Invalid atlas "' segmentation_name '": does not contain any labels.']);
        return;
    end
    
    % Load tissue property
    FilesTissuesProperty = '';
    if isfield(options, 'FilesTissuesProperty')
        FilesTissuesProperty = options.FilesTissuesProperty;
    end
    
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
    
    flag_overwrite_fluences = options.mcxlab_overwrite_fluences;
    
    % Mcxlab configuration
    cfg = struct();
    cfg.gpuid       = options.mcxlab_gpuid;
    cfg.autopilot   = 1;
    cfg.seed        = hex2dec('623F9A9E'); % set seed to make the simulation repeatible
    
    if strcmp(bst_get('OsType'),'mac64arm') 
        % slower but avoid having issue with macos timeout
        cfg.nphoton     = 1e6;
        cfg.respin      = options.mcxlab_nphoton;
    else
        cfg.nphoton     = options.mcxlab_nphoton * 1e6;
        cfg.respin      = 1;  
    end
    cfg.vol         = uint8(sSegmentation.Cube); % segmentation
    cfg.unitinmm    = sSegmentation.Voxsize(1);   % defines the length unit for a grid ( voxel) edge length [1.0]
    cfg.isreflect   = 1; % reflection at exterior boundary
    cfg.isrefint    = 1; % 1-index mismatch at inner boundaries, [0]-matched index
    cfg.tstart      = 0;
    cfg.tend        = 5e-9;
    cfg.tstep       = 5e-9;
    cfg.isnormalized = 1; % [1]-normalize the output flux to unitary source
    cfg.issrcfrom0   = 0; % [0]- first voxel is [1 1 1] ; 1-first voxel is [0 0 0],
    
    tic
    bst_progress('start', 'Compute fluences', sprintf('Computing fluences for %d vertices and %d wavelengths', nb_vertex, nb_wavelengths), 0, nb_vertex * nb_wavelengths);
    
    for iwl = 1:nb_wavelengths
    
        [cfg.prop, errors] = nst_get_tissues_optical_properties(sSegmentation.Labels, wavelengths(iwl), FilesTissuesProperty);
        if ~isempty(errors)
            bst_error(strjoin([ 'Some error occured when extracting the optical properties of the tissues:',  errors], [newline, '---- ']))
            return;
        end
    
        for ivertx =  1:nb_vertex
    
            output_dir = options.outputdir;
            fluence_fn = process_nst_import_head_model('get_fluence_fn', valid_vertices(ivertx), wavelengths(iwl));
    
            if exist(fullfile(output_dir, fluence_fn), 'file') == 2 && flag_overwrite_fluences == 0
                fprintf('Fluence for head vertex %g exists in the target folder\n', valid_vertices(ivertx));
                bst_progress('inc', 1);
    
                continue;
            end
            
            % Initialize position and source direction
            cfg.srcpos = head_vertices_mri(valid_vertices(ivertx), :);    
            cfg.srcdir = head_normals(valid_vertices(ivertx),:);
            
            % Project source position in the volume -- to investigate

            % Original code
            cfg.srcpos = projectPositionInVolume(cfg.vol, cfg.srcpos + 1, cfg.srcdir);
            % Propisition 1: remove this line: keep the original position of the optodes
            % proposition 2: project, but remove the +1; need to add projection on line 346           
            % cfg.srcpos = projectPositionInVolume(cfg.vol, cfg.srcpos + 1, cfg.srcdir);
            % proposition 3: project, move the position following the normal
            % cfg.srcpos = projectPositionInVolume(cfg.vol, cfg.srcpos + cfg.srcdir, cfg.srcdir);

            if strcmp(options.software, 'mcxlab-cuda')
                fluenceRate = mcxlab(cfg); % fluence rate [1/mm2 s]
            else
                fluenceRate = mcxlabcl(cfg); % fluence rate [1/mm2 s]
            end
            
            if isempty(fluenceRate.data) || isinf(fluenceRate.stat.normalizer)
                if strcmp(bst_get('OsType'), 'mac64arm')
                    bst_error('ERROR: MCXlabCL reach time limit (Try again with less photons)');
                else
                    bst_error(sprintf('ERROR:Fluence %d cannot be computed (see command windows)', valid_vertices(ivertx)));
                end
                return;
            end
            
            fluence_vol             = fluenceRate.data .* cfg.tstep;
            fluence_flat_sparse_vol = sparse(double(fluence_vol(:)));

            % Original code
            reference_voxel_index   = cfg.srcpos;
            % If doing proposition 2, we need either:
            % Proposition 2a. Project inside the volume to get the corresponding voxel
            % reference_voxel_index   = projectPositionInVolume(cfg.vol, cfg.srcpos, cfg.srcdir);
            % Or, we do as 3. and we already move along the normal
            % reference_voxel_index   = projectPositionInVolume(cfg.vol, cfg.srcpos + cfg.srcdir , cfg.srcdir);
    
            save(fullfile(output_dir, fluence_fn), 'fluence_flat_sparse_vol','reference_voxel_index');
            bst_progress('inc', 1);
        end
    end
    
    fprintf('\n%d fluence volumes computed in %.2f seconds\n', nb_vertex * nb_wavelengths, toc);
end

function position = projectPositionInVolume(segmentation, position, normals)
%==========================================================================
% Project optodes on non background tissue :
% if the specified optode position is outside the domain, move the optode
% along the initial vector until it hits the domain
%==========================================================================

    % Put the position inside the volume domain (should do nothing)
    if (any(position < 1) || position(1) > size(segmentation, 1) || position(2) > size(segmentation, 2) || position(3) > size(segmentation, 3))
        error('One optodes is outside the volume. Aborting');
    end
    
    % if the voxel coordinates are not integers change them for the nearest integer
    position  = round(position);
    step_size = 1;
    
    %==========================================================================
    % searching non-zero voxel failed along the incident vector
    %==========================================================================
    if segmentation(position(1), position(2), position(3)) > 0  
        % Nothing to do
        return
    end
    
    iter = 0;
    fprintf('Optode at position (%f %f %f) is located outside the domain. ', position);
    
    while segmentation(position(1),position(2),position(3)) == 0 && iter < 400
    
        % Compute the new position and clip
        position = round(position + step_size * normals);
        position = clip(position, [1, 1, 1], size(segmentation));
    
        % Update the position
        iter = iter+1;
    end
    
    fprintf('fixing position to (%f %f %f)\n\n', position);
end

function value = clip(value, min_value, max_value)
% Clip value between min_value and max_value
    for iValue = 1:length(value)
        value(iValue) = max(value(iValue), min_value(iValue));
        value(iValue) = min(value(iValue), max_value(iValue));
    end
end