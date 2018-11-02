function varargout = process_nst_cortical_projection_mne( varargin )

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
% Authors: Thomas Vincent, Alexis Machado (2018)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Cortical projection - MNE';
sProcess.FileTag     = '';
sProcess.Category    = 'File';
sProcess.SubGroup    = 'NIRS - wip';
sProcess.Index       = 1206;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'results', 'results'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.isSeparator = 1;

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};
sStudy = bst_get('Study', sInputs.iStudy);

if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider process "Compute head model from fluence"');
    return;
end
head_model_fn = sStudy.HeadModel(sStudy.iHeadModel).FileName;
head_model = in_bst_headmodel(head_model_fn);

ChanneMat = in_bst_channel(sInputs(1).ChannelFile);

if ~strcmp(head_model.HeadModelType, 'surface')
    bst_error('Projection only works for surface head model');
    return;
end

if ndims(head_model.Gain) ~= 3
    % TODO: better test shape consistency
    bst_error('Bad shape of gain matrix, must be nb_pairs x nb_wavelengths x nb_vertices');
    return;
end

% Load recordings
if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs(1).FileName);
    events = sDataIn.Events;
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
    sDataRaw = in_bst_data(sInputs(1).FileName, 'F');
    events = sDataRaw.F.events;
end

% Separate NIRS channels from others (NIRS_AUX etc.)
[fnirs, fchannel_def, nirs_other, channel_def_other] = ...
    process_nst_mbll('filter_data_by_channel_type', sDataIn.F', ChanneMat, 'NIRS');


% Remove bad channels: they won't enter MBLL computation so no need to keep them
[good_nirs, good_channel_def] = process_nst_mbll('filter_bad_channels', sDataIn.F',...
                                                 ChanneMat, sDataIn.ChannelFlag);
    
[nirs_psig, pair_names, pair_loc, pair_indexes] = process_nst_mbll('group_paired_channels', good_nirs, good_channel_def);
try
    sensitivity_surf = process_nst_import_head_model('get_sensitivity_from_chans', head_model, pair_names);
catch ME  
    if ~isempty(sStudy.iHeadModel)
        bst_warning('Warning', sProcess, sInputs, 'Given head model is not consistent with current montage. Using default one.');
        head_model_fn = sStudy.HeadModel(sStudy.iHeadModel).FileName;
        head_model = in_bst_headmodel(head_model_fn);
        sensitivity_surf = process_nst_import_head_model('get_sensitivity_from_chans', ...
                                                         head_model, pair_names);
    else
        if strcmp(ME.identifier, 'NIRSTORM:HeadmodelMismatch')
            bst_report('Error', sProcess, sInputs, ME.message);
            return;
        else
            rethrow(ME);
        end
    end
end
% nirs_psig (nb_pairs x nb_wavelengths x nb_samples
% nb_nodes = size(sensitivity_surf, 3);
% nb_pairs = size(pair_ichans, 1);

% nb_wavelengths x [HbO, HbR]
ext_coeffs = process_nst_mbll('get_hb_extinctions', ChanneMat.Nirs.Wavelengths);

param.sensors.cov.flag_cov = 1;
param.sensors.cov.window = [sDataIn.Time(1) sDataIn.Time(1)+5];
[pdata_hbo, pdata_hbr] = mfip_inverse_problem(squeeze(nirs_psig(:,1,:)), ...
    squeeze(nirs_psig(:,2,:)), ...
    sDataIn.Time, sensitivity_surf, ...
    ext_coeffs, param);


[sStudy, ResultFile] = nst_bst_add_surf_data(pdata_hbo, sDataIn.Time, ...
    head_model, 'hbo_projected' , ...
    sInputs, sStudy,  ...
    'Projected HbO signals');
OutputFiles{end+1} = ResultFile;

[sStudy, ResultFile] = nst_bst_add_surf_data(pdata_hbr, sDataIn.Time, ...
    head_model, 'hbr_projected' , ...
    sInputs, sStudy,  ...
    'Projected HbR signals');
OutputFiles{end+1} = ResultFile;

clear pdata_hbo pdata_hbr;

%     OutputFiles{end+1} = ResultFile;
% for iw=1:size(sensitivity_surf, 2)
%     mixing_mat = squeeze(sensitivity_surf(:, iw, :));
%     nz = all(mixing_mat ~= 0, 2);
%     mixing_mat(nz,:) = mixing_mat(nz,:) ./ repmat(sum(mixing_mat(nz,:), 2), 1, nb_pairs);
%     projected_fnirs = mixing_mat * sDataIn.F(pair_ichans(:,iw), :);
%     % TODO: comment name including measure tag
%
%     comment_measure = sprintf('wl%d', ChanneMat.Nirs.Wavelengths(iw));
%
%
%     % label = [sDataIn.Comment '_' comment_measure '_projected'];
%     label = [comment_measure '_projected'];
%     [sStudy, ResultFile] = nst_bst_add_surf_data(projected_fnirs, sDataIn.Time, ...
%         head_model, label , ...
%         sInputs.iStudy, sStudy,  ...
%         'Projected fNIRS signals');
%
%     OutputFiles{end+1} = ResultFile;
%
% end

end

function [pdata_hbo, pdata_hbr] = mfip_inverse_problem(data_w1, data_w2, t, mat_A, ex, param)
% data_wX contains signals (nb_pairs x nb_samples)
% mat_A is lead field matrix (nb_channels x nb_wavelengths x nb_voxels)
% ex contains the extinction coeffs nb_wavelength x [HbO HbR]

% A_wX -> (nb_pairs x nb_voxels)
A_w1 = squeeze(mat_A(:,1,:));
A_w2 = squeeze(mat_A(:,2,:));

nSamples = length(t);

nVox = size(mat_A, 3);

% data_w1=data_w1'; % nSensors by nSamples
% data_w2=data_w2';
data_topo = [data_w1 ; data_w2]; % nb_pairs*2 x nSamples

flag_show_l_curve = 0;
param.inverse.Tikkonov.lambda{1} = 0.5;

g=[A_w1*ex(1,1) A_w1*ex(1,2); A_w2*ex(2,1) A_w2*ex(2,2)];
%g=blkdiag(A_w1,A_w2);
I=speye(size(g,1));
S=svd(g); %[U,S,V] = svd(g);
if flag_show_l_curve & numel(param.inverse.Tikkonov.lambda{1}>1)
    % L curve
    %---------------------
    blocks=1:100:nSamples+1; % +1 because we substract one in the next block of code for the upper limit
    if ~(blocks(end)==nSamples+1)
        blocks(end+1)=nSamples+1;
    end
    
    for iL=1:numel(param.inverse.Tikkonov.lambda{1})
        % M=g'/(g*g'+param.inverse.Tikkonov.lambda{1}(iL)*I);
        M=g'/(g*g'+param.inverse.Tikkonov.lambda{1}(iL)*max(S)*I);
        
        J=zeros(nVox*2,nSamples);
        for iB=1:numel(blocks)-1
            fprintf('block %g\n',iB)
            J(:,blocks(iB):blocks(iB+1)-1)= M*data_topo(:,blocks(iB):blocks(iB+1)-1);
        end
        
        norm_res(iL)=norm(data_topo-g*J,2);
        norm_sol(iL)=norm(J,2);
        
    end
    clear M blocks J; % M cleared??
    
    ismonotonic = ismonotonic(norm_res, 1);
    ismonotonic = ismonotonic(norm_sol, 1);
    
    subplot(1,1,1)
    % plot(log10(norm_res),log10(norm_sol),'o');
    plot(norm_res,norm_sol,'o');
    xlabel('||J||2');
    ylabel('||d-gJ||2');
    axis('tight')
    [k_corner,info] = corner(norm_res,norm_sol);
    hold on
    plot(log10(norm_res(k_corner)),log10(norm_sol(k_corner)),'+');
    param.inverse.Tikkonov.lambda{1}(k_corner)
    
else
    M=g'/(g*g'+param.inverse.Tikkonov.lambda{1}(1)*max(S)*I);
end

%TODO: recommand small amount of vertices
%TODO: use sparse storage?
pdata_hb = M * data_topo;
pdata_hbo = pdata_hb(1:nVox, :);
pdata_hbr = pdata_hb((nVox+1):(2*nVox), :);

clear pdata_hb;
end




