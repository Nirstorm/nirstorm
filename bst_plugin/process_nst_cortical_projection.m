function varargout = process_nst_cortical_projection( varargin ) %#ok<STOUT,>

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
sProcess.Comment     = 'Cortical projection';
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

sProcess.options.method.Comment = 'Method';
sProcess.options.method.Type    = 'combobox';
choices = methods();
sProcess.options.method.Value   = {choices.Sensitivity_based_interpolation, ...
                                   fieldnames(choices)};

sProcess.options.compute_hbt.Comment = 'Compute HbT';
sProcess.options.compute_hbt.Type    = 'checkbox';
sProcess.options.compute_hbt.Value   =  0;                        
                               
sProcess.options.sparse_storage.Comment = 'Sparse storage';
sProcess.options.sparse_storage.Type    = 'checkbox';
sProcess.options.sparse_storage.Value   =  0;

sProcess.options.save_mixing_mat.Comment = 'Save mixing matrix';
sProcess.options.save_mixing_mat.Type    = 'checkbox';
sProcess.options.save_mixing_mat.Value   =  0;
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
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
end

% Select valid NIRS channels

montage_info = nst_montage_info_from_bst_channels(ChanneMat.Channel, sDataIn.ChannelFlag);
nirs_sig_wl1 = sDataIn.F(montage_info.pair_ichans(:,1),:);
nirs_sig_wl2 = sDataIn.F(montage_info.pair_ichans(:,2),:);
nb_samples = size(sDataIn.F, 2);

sensitivity_surf = process_nst_import_head_model('get_sensitivity_from_chans', head_model, montage_info.pair_names);
    
% nb_wavelengths x [HbO, HbR]
% cm^-1.l.mol^-1
ext_coeffs = process_nst_mbll('get_hb_extinctions', ChanneMat.Nirs.Wavelengths);

param.sensors.cov.flag_cov = 1;
param.sensors.cov.window = [sDataIn.Time(1) sDataIn.Time(1)+5];

all_methods = methods();
i_method = sProcess.options.method.Value{1};
method_names = fieldnames(all_methods);
mtag = get_method_tag(method_names{i_method});
switch i_method
    case all_methods.Sensitivity_based_interpolation
        [pdata_dOD_proj, mixing_mat] = sensitivity_based_interpolation(nirs_sig_wl1, nirs_sig_wl2, ...
                                                                       sensitivity_surf);
    case all_methods.MNE
        [pdata_dOD_proj, mixing_mat] = mfip_inverse_problem(nirs_sig_wl1, nirs_sig_wl2, ...
                                                            sensitivity_surf, param);        
    otherwise
        error('Invalid method');
end

nb_vertices = size(sensitivity_surf, 3);
hb_types = get_hb_types(sProcess.options.compute_hbt.Value);

pinv_ext = pinv(ext_coeffs);
pdata_hbo = zeros(nb_vertices, nb_samples);
pdata_hbr = zeros(nb_vertices, nb_samples);
for ivox=1:nb_vertices
   hb_sigs = pinv_ext * pdata_dOD_proj([ivox ivox+nb_vertices], :);
   pdata_hbo(ivox,:) = hb_sigs(strcmp(hb_types, 'HbO'), :);
   pdata_hbr(ivox,:) = hb_sigs(strcmp(hb_types, 'HbR'), :);
end

extra.DisplayUnits = 'mol.l-1';
[sStudy, ResultFile] = nst_bst_add_surf_data(pdata_hbo, sDataIn.Time, ...
    head_model, 'hbo_proj', sprintf('HbO cortex (%s)', mtag), ...
    sInputs, sStudy,  'Projected HbO signals', [], sProcess.options.sparse_storage.Value, ...
    extra);
OutputFiles{end+1} = ResultFile;

[sStudy, ResultFile] = nst_bst_add_surf_data(pdata_hbr, sDataIn.Time, ...
    head_model, 'hbr_proj', sprintf('HbR cortex (%s)', mtag) , ...
    sInputs, sStudy, 'Projected HbR signals', [], sProcess.options.sparse_storage.Value, ...
    extra);
OutputFiles{end+1} = ResultFile;

if sProcess.options.compute_hbt.Value
    pdata_hbt = pdata_hbo + pdata_hbr;
    [sStudy, ResultFile] = nst_bst_add_surf_data(pdata_hbt, sDataIn.Time, ...
                                                 head_model, 'hbt_proj', sprintf('HbT cortex (%s)', mtag), ...
                                                 sInputs, sStudy,  'Projected HbT signals', [], sProcess.options.sparse_storage.Value, ...
                                                 extra);
    OutputFiles{end+1} = ResultFile;
    clear pdata_hbt;
end

clear pdata_hbo pdata_hbr;

if sProcess.options.save_mixing_mat.Value
    sCortex = in_tess_bst(head_model.SurfaceFile);
    nb_vertices = size(sCortex.Vertices, 1);
    % TODO: remap treated channels into original channels listing
    i_chans = [montage_info.pair_ichans(:,1) ; montage_info.pair_ichans(:,2)];
    chan_data = zeros(size(sDataIn.F,1), nb_vertices);
    chan_data(i_chans, :) = mixing_mat';
    save_chan_data(chan_data, sDataIn, sInputs.iStudy, sStudy, sprintf('Mixings (%s)', mtag), 'mixing');
end

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

function hb_types = get_hb_types(enable_hbt)
if nargin < 1
    enable_hbt = 0;
end
% Return Hb type order as used for the cols of the matrix of exctinction coeffs
if enable_hbt
    hb_types = {'HbO', 'HbR', 'HbT'};
else
    hb_types = {'HbO', 'HbR'};
end
end

function save_chan_data(chan_data, sDataIn, iStudy, sStudy, comment, tag)
    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = chan_data;
    sDataOut.Comment      = comment;
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = 1:size(chan_data,2);
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = 1;
    sDataOut.Events       = [];
    sDataOut.DisplayUnits = 'U.A.';

    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ['data_' tag]);
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudy, OutputFile, sDataOut);
end

function [pdata_dOD_proj, mixings] = mfip_inverse_problem(data_w1, data_w2, mat_A, param)
% data_wX contains signals (nb_pairs x nb_samples)
% mat_A is lead field matrix (nb_channels x nb_wavelengths x nb_vertices)
% ex contains the extinction coeffs nb_wavelength x [HbO HbR]
%
% Outputs:
%     - pdata_dOD_proj, size (nb_vertices*nb_wls) x nb_samples
%     - M: mixing matrix, size: (nb_vertices*2) x (nb_pairs*nb_wls)

% A_wX -> (nb_pairs x nb_vertices)
[nb_pairs, nb_wls, nb_vertices] = size(mat_A); %#ok<ASGLU>
A_w1 = zeros(nb_pairs, nb_vertices);
A_w2 = zeros(nb_pairs, nb_vertices);
A_w1(:,:) = mat_A(:,1,:);
A_w2(:,:) = mat_A(:,2,:);

% data_w1=data_w1'; % nSensors by nSamples
% data_w2=data_w2';
data_topo = [data_w1 ; data_w2]; % (nb_pairs*nb_wls) x nSamples

flag_show_l_curve = 0;
param.inverse.Tikkonov.lambda{1} = 0.1;

%g=[A_w1*ex(1,1) A_w1*ex(1,2); A_w2*ex(2,1) A_w2*ex(2,2)]; % (nb_pairs*nb_wls) x (nb_vertices*2)

g=blkdiag(A_w1,A_w2); % (nb_pairs*nb_wls) x (nb_vertices*nb_wls)
I=speye(size(g,1));
% S=svd(g); %[U,S,V] = svd(g);
S = svd(g*g');
if flag_show_l_curve && numel(param.inverse.Tikkonov.lambda{1}) > 1
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
    % (nb_vertices*2) x (nb_pairs*nb_wls) 
    M=g'/(g*g' + param.inverse.Tikkonov.lambda{1}(1)*max(S)*I);
end

%TODO: recommand small amount of vertices 
pdata_dOD_proj = M * data_topo; % (nb_vertices*nb_wls) x nb_samples

mixings = M(1:nb_vertices, :);

end

function [pdata_dOD_proj, M] = sensitivity_based_interpolation(data_w1, data_w2, mat_A)
% data_wX contains signals (nb_pairs x nb_samples)
% mat_A is lead field matrix (nb_channels x nb_wavelengths x nb_vertices)

nb_samples = size(data_w1, 2);
data_topo = {data_w1, data_w2};

[nb_pairs, nb_wls, nb_vertices] = size(mat_A);

A_wl = zeros(nb_pairs, nb_vertices);
pdata_dOD_proj = zeros(nb_vertices*nb_wls, nb_samples);
M = zeros(nb_vertices, nb_pairs*nb_wls);
for iwl=1:nb_wls
    A_wl(:,:) = mat_A(:,iwl,:);
    M_wl = (A_wl ./ repmat(sum(A_wl,2), 1, nb_vertices))';
    M_wl(isnan(M_wl)) = 0;
    pdata_dOD_proj((nb_vertices*(iwl-1)+1):(nb_vertices*iwl), :) = M_wl * data_topo{iwl};
    M(:, (nb_pairs*(iwl-1)+1):(nb_pairs*iwl)) = M_wl;
end

end

function enum = methods()
enum.Sensitivity_based_interpolation = 1;
enum.MNE = 2;
end

function tag = get_method_tag(method_name)

if length(method_name) <= 3
    tag = method_name;
else
    tag = strjoin(cellfun(@(s) s(1), strsplit('Sensitivity_based_interpolation', '_'), ...
                    'UniformOutput', 0), ...
                  '');
end

end


