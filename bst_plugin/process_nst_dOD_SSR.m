function varargout = process_nst_dOD_SSR( varargin )

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
% Authors: Thomas Vincent (2015-2016)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'MBLL - OD to delta OD_SSR';
    sProcess.FileTag     = ' | dOD';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1004; %0: not shown, >0: defines place in the list of processes
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSEEGVisualCheckerboard';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'}; 
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options

    sProcess.options.option_baseline_method.Comment = 'Baseline method';
    sProcess.options.option_baseline_method.Type    = 'combobox';
    sProcess.options.option_baseline_method.Value   = {1, {'mean', 'median'}}; % {Default index, {list of entries}}
    
    sProcess.options.option_do_SuperficalRegression.Comment = 'Superfical Regression';
    sProcess.options.option_do_SuperficalRegression.Type    = 'checkbox';
    sProcess.options.option_do_SuperficalRegression.Value   = 0;
    
    sProcess.options.option_Superfical_Channel.Comment = 'Superfical Channel [coma-separated list]';
    sProcess.options.option_Superfical_Channel.Type    = 'text';
    sProcess.options.option_Superfical_Channel.Value   = ''; 
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>

    % Get option values
    blm_idx         = sProcess.options.option_baseline_method.Value{1};
    baseline_method = sProcess.options.option_baseline_method.Value{2}{blm_idx};
    
    % Load channel file
    ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
    do_SuperficalRegression = sProcess.options.option_do_SuperficalRegression.Value;
    
    
    %load superfical channel info
    if do_SuperficalRegression
    superfical_chann = sProcess.options.option_Superfical_Channel.Value;
    try
        supf_chann = textscan(superfical_chann, '%s','Delimiter',',');
    catch
        bst_report('Error', sProcess, [], 'Superfical Channels must be integers separated by comas');
        return
    end
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
    
    % Remove bad channels: they won't enter dOD computation so no need to keep them     
    % Separate NIRS channels from others (NIRS_AUX etc.)
    to_keep = sDataIn.ChannelFlag ~= -1 & strcmpi({ChanneMat.Channel.Type}, 'NIRS')';
    
    
    % Remove bad channels: they won't enter MBLL computation so no need to keep them 
    [good_nirs, good_channel_def] = process_nst_mbll_SSR('filter_bad_channels',sDataIn.F', ChanneMat, sDataIn.ChannelFlag);
    
    % Separate NIRS channels from others (NIRS_AUX etc.)                                                
    [fnirs, fchannel_def, nirs_other, channel_def_other] = ...
        process_nst_mbll_SSR('filter_data_by_channel_type',good_nirs, good_channel_def, 'NIRS');
    
    
    % Apply dOD computation
    if do_SuperficalRegression
    [nirs_dOD, channel_dOD_def] = Compute(fnirs, baseline_method, fchannel_def, do_SuperficalRegression, supf_chann);
    else
    [nirs_dOD, channel_dOD_def] = Compute(fnirs, baseline_method, fchannel_def, do_SuperficalRegression);
    end
    
    [final_nirs, ChannelMat] = process_nst_mbll_SSR('concatenate_data',nirs_dOD, channel_dOD_def, nirs_other, channel_def_other);
    
    % Create new condition because channel definition is different from original one
     cond_name = sInputs.Condition;
     if strcmp(cond_name(1:4), '@raw')
        cond_name = cond_name(5:end);
     end
     if do_SuperficalRegression
         iStudy = db_add_condition(sInputs.SubjectName, [cond_name, '_dOD_SSR']);
     else
         iStudy = db_add_condition(sInputs.SubjectName, [cond_name, '_dOD']);
     end
     sStudy = bst_get('Study', iStudy);
%     
%     % Save channel definition
     [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
     db_set_channel(iChannelStudy, ChannelMat, 0, 0);
        
    % Save time-series data
%     final_nirs = sDataIn.F;
%     final_nirs(to_keep, :) = final_nirs';
    sDataOut = db_template('data');
    sDataOut.F            = final_nirs';
    sDataOut.Comment      = 'NIRS dOD';
    %sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.ChannelFlag  = ones(size(final_nirs, 2), 1);
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = 1;
    sDataOut.Events       = events;
    sDataOut.DisplayUnits = 'delta OD';

    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_hb');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudy, OutputFile, sDataOut);
end


function [nirs_dOD, channel_dOD_def] = Compute(nirs_sig, normalize_method, channel_def, do_SuperficalRegression, supf_chann)
%% Normalize given nirs signal
% Args:
%    - nirs_sig: matrix of double, size:  nb_samples x nb_channels
%        NIRS signal to normalize
%   [- method]: str, choices are: 'mean' and 'median', default is 'mean'
%        Normalization method.
%        * 'mean': divide given nirs signal by its mean, for each wavelength
%        * 'median': divide given nirs signal by its median, for each wavelength
% 
% Output: matrix of double, size: nb_samples x nb_channels
%    Normalized NIRS signal.

switch normalize_method
    case 'mean'
        od_ref = mean(nirs_sig, 1);
    case 'median'
        od_ref = median(nirs_sig, 1);
end

[nirs_psig, pair_names, pair_loc, pair_indexes] = process_nst_mbll_SSR('group_paired_channels',nirs_sig, channel_def);
pair_distances = process_nst_mbll_SSR('cpt_distances',channel_def.Channel, pair_indexes) .* 100; %convert to cm
nb_samples = size(nirs_sig, 1);
nb_pairs = length(pair_names);
dOD = zeros(nb_pairs, 2, nb_samples);

if do_SuperficalRegression
    [~,idx_supf_chann] = find(ismember(pair_names,supf_chann{1}'));
    delta_od_supf = zeros(2,nb_samples);
    for i_supf_chann = 1:length(idx_supf_chann)
       delta_od_supf = process_nst_mbll_SSR('normalize_nirs',squeeze(nirs_psig(idx_supf_chann(i_supf_chann), :, :)), ...
                              normalize_method)+delta_od_supf; 
       disp(['Superfical Channel-' pair_names(idx_supf_chann(i_supf_chann)) num2str(pair_distances(idx_supf_chann(i_supf_chann)),'%.2f') 'cm']);
    end
    delta_od_supf = delta_od_supf./length(idx_supf_chann);
end

for ipair=1:size(nirs_psig, 1)
    
    delta_od = process_nst_mbll_SSR('normalize_nirs',squeeze(nirs_psig(ipair, :, :)), ...
                              normalize_method);
    
    % SSR regression  Gregg 2010 Frontiers
    if do_SuperficalRegression
        a_data = dot(delta_od',delta_od_supf')./dot(delta_od_supf',delta_od_supf');
        delta_od = delta_od - repmat(a_data,[size(delta_od_supf,2),1])'.*delta_od_supf;
        dOD(ipair, 1:2, :) = delta_od;
    else
        dOD(ipair, 1:2, :) = delta_od;
    end
    
end
[nirs_dOD, channel_dOD_def] = pack_dOD_channels(dOD, pair_names, pair_loc, channel_def);
%nb_samples = size(nirs_sig, 1);

%delta_od = -log( nirs_sig ./ repmat(od_ref, nb_samples, 1) );
end

% function [fdata, fchannel_def] = ...
%     filter_bad_channels(data, channel_def, channel_flags)
% %% Filter the given data based on bad channel flags
% %
% % Args:
% %    - data: matrix of double, size: nb_samples x nb_channels
% %      data time-series to filter
% %    - channel_def: struct
% %        Defintion of channels as given by brainstorm
% %        Used field: Channel
% %    - channel_flags: array of int, default: []
% %        Channel flags. Channel with flag -1 are filtered
% %        If [] is given, all channels are kept.
% if nargin < 3 || isempty(channel_flags)
%     channel_flags = ones(length(channel_def.Channel), 1);
% end
% 
% kept_ichans = channel_flags' ~= -1;
% 
% fchannel_def = channel_def;
% fchannel_def.Channel = channel_def.Channel(kept_ichans);
% fdata = data(:, kept_ichans);
% end

function [nirs_dOD, channel_hb_def] = pack_dOD_channels(nirs, pair_names, pair_loc, channel_def_orig)
%% Reshape given nirs data to have size: time x nb_channels
%% Build new dOD-specific channel listing
%
% Args:
%    - nirs: array of double, size: nb_pairs x 3 x nb_samples
%        NIRS data to reshape (second axis corresponds to [HbO, HbR, HbT])
%    - pair_names: cell array of str, size; nb_pairs
%        Pair names (format: SXDY, where X and Y are the src and det indexes, resp.)
%    - pair_loc: matrix of double, size: (nb_pairs, 3)
%        Spatial coordinates of pairs. Loc(:,1) is source, Loc(:,2) is detector
%    - channel_def_orig: struct
%        Initial defintion of channels as given by brainstorm
%        Used to import field that remain the same after MBBL
%        -> only fields 'Nirs' and 'Channel' are redefined
%
% Outputs:
%   - nirs_hb: matrix of double, size: time x (nb_channels/nb_wavelengths)*2
%       HbO and HbO delta concentration time-series
%   - channel_hb_def: struct 
%       definition of new Hb-related channels. Relevant content:
%           Nirs.Hb = {'HbO', 'HbR', 'HbT'};
%           Channel(ichan1).Name = 'SXDXHbO'; % pair SXDX, HbO component
%           Channel(ichan1).Loc = [C C C]; % pair localization (imported from
%                                          % input channel_def)
%           Channel(ichan1).Group = 'HbO'; 
%           Channel(ichan2).Name = 'SXDXHbR'; % pair SXDX, HbO component
%           Channel(ichan2).Loc = [C C C]; % pair localization (imported from
%                                          % input channel_def, same as paired 
%                                          % channel ichan1)
%           Channel(ichan2).Group = 'HbR'; 
%           Channel(ichan3).Name = 'SXDXHbT'; % pair SXDX, HbT component
%           Channel(ichan3).Loc = [C C C]; % pair localization (imported from
%                                          % input channel_def, same as paired 
%                                          % channel ichan1)
%           Channel(ichan3).Group = 'HbT'; 

channel_hb_def = channel_def_orig;
ichan = 1;
%hb_names = {'HbO', 'HbR'};
for ipair=1:size(nirs, 1)
    for idOD=1:length(channel_hb_def.Nirs.Wavelengths)
        nirs_dOD(:, ichan) = squeeze(nirs(ipair, idOD, :));
         Channel(ichan).Name = [pair_names{ipair} 'WL' num2str(channel_hb_def.Nirs.Wavelengths(idOD))];
         Channel(ichan).Loc = squeeze(pair_loc(ipair, :, :));
         Channel(ichan).Group = ['WL' num2str(channel_hb_def.Nirs.Wavelengths(idOD))];
         Channel(ichan).Comment = [];
         Channel(ichan).Orient = [];
         Channel(ichan).Weight = 1;
         Channel(ichan).Type = 'NIRS';
        ichan = ichan + 1;
    end
end

%channel_hb_def.Nirs = rmfield(channel_hb_def.Nirs, 'Wavelengths');

%channel_hb_def.Nirs.Hb = hb_names;
channel_hb_def.Channel = Channel;
channel_hb_def.Comment = ['NIRS-BRS sensors (' num2str(length(Channel)) ')'];
end

