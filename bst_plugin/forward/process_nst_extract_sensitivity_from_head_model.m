function varargout = process_nst_extract_sensitivity_from_head_model( varargin )

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
% Authors: Edouard Delaire (2023) 
%          Thomas Vincent, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract sensitivity surfaces from head model';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1405;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;


    % === Process description
    sProcess.options.label1.Comment = ['<b>Light sensitivity scale.</b> <BR>' ...
                                        'Export the nirs global and channel specific sensitivity based on the head model<BR>', ...
                                       'Note: if the sensitivity is exported in the log scale (db), the map is thresholded at -2db.'];
    sProcess.options.label1.Type = 'label';


    sProcess.options.method.Comment = {['Linear'], ...
                                       ['Scale with the global max (dB): <FONT color=#7F7F7F>&nbsp;&nbsp;&nbsp;' ...
                                            'sensitivity = log10(x / max(sensitivity))'],... 
                                       ['Scale with the local max (dB) (only for channel specific sensitivity map): <FONT color=#7F7F7F>&nbsp;&nbsp;&nbsp;' ...
                                            'For each channel, sensitivity = log10(x / max(sensitivity))</FONT>'] ;...
                                       'linear', 'db_global','db_local'};
    sProcess.options.method.Type    = 'radio_label';
    sProcess.options.method.Value   = 'linear';


    sProcess.options.label2.Comment = ['<b>Measure of overlap</b> <BR>'];
    sProcess.options.label2.Type = 'label';

    sProcess.options.export_overlap.Comment = 'Export measure of overlap';
    sProcess.options.export_overlap.Type    = 'checkbox';
    sProcess.options.export_overlap.Value   = 1;    


    
    sProcess.options.label3.Comment = ['<b>NIRS field of view (FOV)</b> <BR>' ...
                                       'Export the FOV used for source reconstruction as new atlas'];
    sProcess.options.label3.Type = 'label';

    sProcess.options.export_FOV.Comment = 'Export NIRS FOV';
    sProcess.options.export_FOV.Type    = 'checkbox';
    sProcess.options.export_FOV.Value   = 1;    
    sProcess.options.export_FOV.Controller = 'FOV';

    % Definition of the options
    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {3, 'cm',2};
    sProcess.options.thresh_dis2cortex.Class   = 'FOV';


end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

% Save the new head model
sStudy = bst_get('Study', sInputs.iStudy);

if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider process "Compute head model from fluence"');
    return;
end

bst_chan_data   = load(file_fullpath(sInputs.FileName), 'ChannelFlag');
ChannelFlag     = bst_chan_data.ChannelFlag;

head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);
sCortex = in_tess_bst(head_model.SurfaceFile);



if ~strcmp(head_model.HeadModelType, 'surface')
    bst_error('Extraction only works for surface head model');
    return;
end


if ndims(head_model.Gain) ~= 3
   bst_error('Bad shape of gain matrix, must be nb_pairs x nb_wavelengths x nb_vertices');
   return;
end

sSubject    = bst_get('Subject', sInputs.SubjectName);
voronoi_fn  = process_nst_compute_voronoi('get_voronoi_fn', sSubject);
    
if ~exist(voronoi_fn, 'file')
    error('Could not find the required Voronoi file.');
end

 
montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel,ChannelFlag);

src_coords = montage_info.src_pos;
det_coords = montage_info.det_pos;

nb_sources          = size(src_coords, 1);
nb_dets             = size(det_coords, 1);

nb_nodes            = size(head_model.Gain   , 3);
nb_Wavelengths      = size(head_model.Gain   , 2);

time        = 1:(nb_sources*100 + nb_dets + 1);
isUsedTime  = zeros(1, length(time));

sensitivity_surf        = zeros(nb_nodes, nb_Wavelengths, length(time));
sensitivity_surf_sum    = zeros(nb_nodes, nb_Wavelengths);

%% Compute sensitivity
for iwl=1:nb_Wavelengths
    swl = ['WL' num2str(ChannelMat.Nirs.Wavelengths(iwl))];
    selected_chans = strcmpi({ChannelMat.Channel.Group}, swl) & (ChannelFlag>0)';
    idx_chan       = find(selected_chans);
        
    sensitivity     = nst_headmodel_get_gains(head_model,iwl, ChannelMat.Channel,idx_chan );


    for iChan = 1:length(idx_chan)
        chan = ChannelMat.Channel(idx_chan(iChan));
        [src_id, det_id] = nst_unformat_channel(chan.Name );

        sensitivity_surf(:,iwl, det_id + src_id*100) = squeeze(sensitivity(iChan,:));
        isUsedTime(det_id + src_id*100)              = 1;
    end

    sensitivity_surf_sum(:,iwl) = sum(sensitivity,  1) ;
end


%% Normalize values and threshold 
threshold_value = -2; % in db

if contains(sProcess.options.method.Value,'db')
    for iwl=1:nb_Wavelengths
        sensitivity_surf_sum(:,iwl) = log10(sensitivity_surf_sum(:,iwl) ./ ( eps + max(sensitivity_surf_sum(:,iwl))));
        mask = zeros(size(sensitivity_surf_sum));
        mask(:,iwl) = sensitivity_surf_sum(:,iwl) < threshold_value;
        sensitivity_surf_sum(mask == 1) = 0;

        k = zeros(1,  size(sensitivity_surf,3));

        if strcmp(sProcess.options.method.Value,'db_local') % channel wise 
            k(1,:) =  squeeze(max(sensitivity_surf(:, iwl, :))) + eps;     
        else % Global normalisation 
            k(1,:) = max(max(sensitivity_surf(:, iwl, :))) + eps;   
        end    
        
        sensitivity_surf(:,iwl,isUsedTime == 1) = log10( squeeze(sensitivity_surf(:, iwl, isUsedTime == 1)) ./ repmat(k(isUsedTime == 1),nb_nodes,1));
        
        mask = zeros(size(sensitivity_surf));
        mask(:,iwl,:) = sensitivity_surf(:,iwl,:) < threshold_value;
        sensitivity_surf(mask == 1 ) = 0;
    end
else
    for iwl=1:nb_Wavelengths
        mask = zeros(size(sensitivity_surf_sum));
        mask(:,iwl) = sensitivity_surf_sum(:,iwl) < 10^(threshold_value)*max(sensitivity_surf_sum(:,iwl));
        sensitivity_surf_sum(mask == 1) = 0;

        mask = zeros(size(sensitivity_surf));
        mask(:,iwl,:) = sensitivity_surf(:,iwl,:) <  10^(threshold_value)*max(sensitivity_surf(:,iwl,:),[],'all');
        sensitivity_surf(mask == 1 ) = 0;
    end
end   

%% Save sensitivity 
for iwl=1:size(sensitivity_surf, 2)


    [sStudy, ResultFile] = add_surf_data(repmat(squeeze(sensitivity_surf_sum(:,iwl)), [1,2]), [0 1], ...
                                         head_model, ['Summed sensitivities - WL' num2str(iwl)], ...
                                         sInputs.iStudy, sStudy,  ...
                                         'sensitivity imported from MCXlab');
        
    OutputFiles{end+1} = ResultFile;


    [sStudy, ResultFile] = add_surf_data( squeeze(sensitivity_surf(:,iwl,:)), time, ...
        head_model, ['Sensitivities - WL' num2str(iwl)], ...
        sInputs.iStudy, sStudy, 'sensitivity imported from MCXlab');
    OutputFiles{end+1} = ResultFile;
end
%% Estimate Coverage

%threshold for coverage
p_thresh = 1;
act_vol = 1000; % A definir comme un parametre donne par l'utilisateur
        
sVoronoi = in_mri_bst(voronoi_fn);

median_voronoi_volume = process_nst_compute_voronoi('get_median_voronoi_volume', sVoronoi);  
delta_mu_a = 0.001;
threshold = compute_threshold(p_thresh, act_vol, median_voronoi_volume, delta_mu_a);

coverage = sensitivity_surf > threshold ;

for iwl=1:size(sensitivity_surf, 2)

    [sStudy, ResultFile] = add_surf_data(repmat(squeeze(sum(coverage(:, iwl,:), 3)), [1,2]), [0 1], ...
                                         head_model, ['Summed Coverage - WL' num2str(iwl)], ...
                                         sInputs.iStudy, sStudy,  ...
                                         'sensitivity imported from MCXlab');
        
    OutputFiles{end+1} = ResultFile;

    [sStudy, ResultFile] = add_surf_data( ...
                                         squeeze(coverage(:,iwl,:)), ...
                                         time, ...
                                        head_model, ...
                                        ['Coverage - WL' num2str(iwl)], ...
                                        sInputs.iStudy, ...
                                        sStudy, ...
                                        'sensitivity imported from MCXlab');
    OutputFiles{end+1} = ResultFile;
end


%% Estimate overlaps - Display sensitivity as function of number of overlap

thresholds = linspace(min(sensitivity_surf,[],'all'), max(sensitivity_surf,[],'all'), 10);


for iwl=1:size(sensitivity_surf, 2)
    tmp = squeeze(sensitivity_surf(:,iwl,isUsedTime == 1));

    overlaps = zeros(size(tmp));
    for iVertex = 1:size(tmp, 1)
        a  = sort(tmp(iVertex,:)  ,2,'descend');

        if contains(sProcess.options.method.Value,'db')
            overlaps(iVertex,1:length(find(a<0)))  = a(a < 0); 
        else
            overlaps(iVertex,1:length(find(a>0)))  = a(a > 0); 
        end
    end
    
    overlaps_number = zeros(size(tmp,1), length(thresholds));
    for iThreshold = 1:length(thresholds)
        if contains(sProcess.options.method.Value,'db')
            overlaps_number(:, iThreshold) = sum( tmp > thresholds(iThreshold)  & tmp < 0  ,2 );
        else
            overlaps_number(:, iThreshold) = sum( tmp > thresholds(iThreshold)  & tmp > 0  ,2 );
        end
    end

    if sProcess.options.export_overlap.Value 

        [sStudy, ResultFile] = add_surf_data(overlaps, 1:size(tmp,2), ...
                                             head_model, sprintf('Overlap WL %d (sensitivity)',iwl), ...
                                             sInputs.iStudy, sStudy,  ...
                                             'sensitivity imported from MCXlab');
        OutputFiles{end+1} = ResultFile;

        [sStudy, ResultFile] = add_surf_data(overlaps_number, thresholds, ...
                                         head_model, sprintf('Overlap WL %d (#overlaps)',iwl), ...
                                         sInputs.iStudy, sStudy,  ...
                                         'sensitivity imported from MCXlab');
        OutputFiles{end+1} = ResultFile;

    end
end


%% define the reconstruction FOV
thresh_dis2cortex           = sProcess.options.thresh_dis2cortex.Value{1}*0.01;
[valid_nodes,dis2cortex]    = nst_headmodel_get_FOV(ChannelMat, sCortex, thresh_dis2cortex, ChannelFlag);

if any(strcmp({sCortex.Atlas.Name},'NIRS-FOV'))
    iAtlas = find(strcmp({sCortex.Atlas.Name},'NIRS-FOV'));
else
    sCortex.Atlas(end+1).Name = 'NIRS-FOV';
    iAtlas = length( sCortex.Atlas);
end


sCortex.Atlas(iAtlas).Scouts(end+1)              = db_template('Scout'); 
sCortex.Atlas(iAtlas).Scouts(end).Vertices       = valid_nodes;
sCortex.Atlas(iAtlas).Scouts(end).Seed           = valid_nodes(1);
sCortex.Atlas(iAtlas).Scouts(end).Label          = sprintf('NIRS FOV (%d cm)',sProcess.options.thresh_dis2cortex.Value{1} );
sCortex.Atlas(iAtlas).Scouts(end)                = panel_scout('SetColorAuto',sCortex.Atlas(iAtlas).Scouts(end), length(sCortex.Atlas(iAtlas).Scouts));

bst_save(file_fullpath(head_model.SurfaceFile), sCortex)

% [sStudy, ResultFile] = add_surf_data(repmat(dis2cortex*100, [1,2]), [0 1], ...
%                                  head_model, 'Distance to cortex', ...
%                                  sInputs.iStudy, sStudy,  ...
%                                  'sensitivity imported from MCXlab');


end

function threshold = compute_threshold(p_thresh, act_vol, V_hat, delta_mu_a)
% @========================================================================
% compute_threshold computes the threshold to determine the coverage matrix 
% Formula used : threshold = log((100 + p_thresh) / 100) / ((act_vol / V_hat) * delta_mu_a)
% ========================================================================@
        
        numerator = log10((100 + p_thresh) / 100);
        denomimator = (act_vol / V_hat) * delta_mu_a;
        threshold = numerator / denomimator;
end

function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
                                              iStudy, sStudy, history_comment)
                                          
%% Save a cortical map to brainstorm with given data

    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                            ['results_NIRS_' nst_protect_fn_str(name)]);
                        
    % ===== CREATE FILE STRUCTURE =====
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    ResultsMat.ImageGridAmp = data;
    ResultsMat.Time          = time;
    ResultsMat.DataFile      = [];
    ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName; 
    ResultsMat.HeadModelType = head_model.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = file_short(head_model.SurfaceFile);
    ResultsMat.GridLoc    = [];
    ResultsMat.GridOrient = [];
    ResultsMat.nAvg      = 1;
    % History
    ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
    % Save new file structure
    bst_save(ResultFile, ResultsMat, 'v6');
    % ===== REGISTER NEW FILE =====
    % Create new results structure
    newResult = db_template('results');
    newResult.Comment       = name;
    newResult.FileName      = file_short(ResultFile);
    newResult.DataFile      = ''; %sInputs.FileName;
    newResult.isLink        = 0;
    newResult.HeadModelType = ResultsMat.HeadModelType;
    % Add new entry to the database
    iResult = length(sStudy.Result) + 1;
    sStudy.Result(iResult) = newResult;
    % Update Brainstorm database
    bst_set('Study', iStudy, sStudy);
end
