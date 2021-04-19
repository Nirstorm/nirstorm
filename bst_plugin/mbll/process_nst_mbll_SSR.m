function varargout = process_nst_mbll_SSR( varargin )

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
    sProcess.Comment     = 'MBLL - raw to delta [HbO], [HbR] & [HbT]';
    sProcess.FileTag     = '_Hb';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Work in progress'};
    sProcess.Index       = 1810; %0: not shown, >0: defines place in the list of processes
    sProcess.isSeparator = 1;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Compute_.5BHb.5D_variations_-_Modified_Beer-Lambert_Law';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'}; %TODO: 'raw' -> 'raw' or 'raw' -> 'data'?
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    sProcess.options.option_age.Comment = 'Age';
    sProcess.options.option_age.Type    = 'value';
    sProcess.options.option_age.Value   = {25, 'years', 2};
    
    sProcess.options.option_pvf.Comment = 'PVF';
    sProcess.options.option_pvf.Type    = 'value';
    sProcess.options.option_pvf.Value   = {50, '', 0};
    
    sProcess.options.option_baseline_method.Comment = 'Baseline method';
    sProcess.options.option_baseline_method.Type    = 'combobox';
    sProcess.options.option_baseline_method.Value   = {1, {'mean', 'median'}};    % {Default index, {list of entries}}

    sProcess.options.option_do_plp_corr.Comment = 'Light path length correction';
    sProcess.options.option_do_plp_corr.Type    = 'checkbox';
    sProcess.options.option_do_plp_corr.Value   = 1;
    
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
       
    % TODO: check for negative values
    % Get option values   
    age             = sProcess.options.option_age.Value{1};
    blm_idx         = sProcess.options.option_baseline_method.Value{1};
    baseline_method = sProcess.options.option_baseline_method.Value{2}{blm_idx};
    do_plp_corr     = sProcess.options.option_do_plp_corr.Value;
    do_SuperficalRegression = sProcess.options.option_do_SuperficalRegression.Value;
    pvf = sProcess.options.option_pvf.Value{1};
    
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
    
    % Load channel file
    ChanneMat = in_bst_channel(sInputs(1).ChannelFile);
    %% TODO : ChanneMat.Nirs should be 1xN instead of Nx1 should be fixed when loading data
    if size(ChanneMat.Nirs.Wavelengths,1)>size(ChanneMat.Nirs.Wavelengths,2)
        ChanneMat.Nirs.Wavelengths = ChanneMat.Nirs.Wavelengths';
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
    
    if any(any(sDataIn.F(sDataIn.ChannelFlag~=-1, :) < 0))
        msg = 'Good channels contains negative values. Consider running NISTORM -> Set bad channels';
        bst_error(msg, '[Hb] quantification', 0);
        return;
    end
    
    % Remove bad channels: they won't enter MBLL computation so no need to keep them 
    [good_nirs, good_channel_def] = filter_bad_channels(sDataIn.F', ChanneMat, sDataIn.ChannelFlag);
    
    % Separate NIRS channels from others (NIRS_AUX etc.)                                                
    [fnirs, fchannel_def, nirs_other, channel_def_other] = ...
        filter_data_by_channel_type(good_nirs, good_channel_def, 'NIRS');
    
    % Apply MBLL
    if do_SuperficalRegression
    [nirs_hb, channels_hb] = Compute(fnirs, fchannel_def, age, baseline_method, do_plp_corr, pvf, do_SuperficalRegression, supf_chann); 
     else
     [nirs_hb, channels_hb] = Compute(fnirs, fchannel_def, age, baseline_method, do_plp_corr, pvf, do_SuperficalRegression);    
    end
    % Re-add other channels that were not changed during MBLL
    [final_nirs, ChannelMat] = concatenate_data(nirs_hb, channels_hb, nirs_other, channel_def_other);

    % Create new condition because channel definition is different from original one
    cond_name = sInputs.Condition;
    if strcmp(cond_name(1:4), '@raw')
        cond_name = cond_name(5:end);
    end
    iStudy = db_add_condition(sInputs.SubjectName, [cond_name, '_Hb']);
    sStudy = bst_get('Study', iStudy);
    
    % Save channel definition
    [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudy);
    db_set_channel(iChannelStudy, ChannelMat, 0, 0);
    
    % Save time-series data
    sDataOut = db_template('data');
    sDataOut.F            = final_nirs'; % TOCHECK brainstorm expects MILLIMOL!!
    sDataOut.Comment      = 'Hb [Topo]';
    sDataOut.ChannelFlag  = ones(size(final_nirs, 2), 1);
    sDataOut.Time         = sDataIn.Time;
    sDataOut.DataType     = 'recordings'; 
    sDataOut.nAvg         = 1;
    sDataOut.Events       = events;
    sDataOut.DisplayUnits = 'mol.l-1';

    % Generate a new file name in the same folder
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_hb');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(iStudy, OutputFile, sDataOut);
end

function [fdata, fchannel_def] = ...
    filter_bad_channels(data, channel_def, channel_flags)
%% Filter the given data based on bad channel flags
%
% Args:
%    - data: matrix of double, size: nb_samples x nb_channels
%      data time-series to filter
%    - channel_def: struct
%        Defintion of channels as given by brainstorm
%        Used field: Channel
%    - channel_flags: array of int, default: []
%        Channel flags. Channel with flag -1 are filtered
%        If [] is given, all channels are kept.
if nargin < 3 || isempty(channel_flags)
    channel_flags = ones(length(channel_def.Channel), 1);
end

kept_ichans = channel_flags' ~= -1;

fchannel_def = channel_def;
fchannel_def.Channel = channel_def.Channel(kept_ichans);
fdata = data(:, kept_ichans);
end


function [fdata, fchannel_def, data_other, channel_def_other] = ...
    filter_data_by_channel_type(data, channel_def, channel_types)

%    - channel_types: str or cell array of str
%        Channel types to keep.
%        If [] is given then all channels are kept.

if nargin < 3 || isempty(channel_types)
    channel_types = unique({channel_def.Channel.Type});
else
    if isstr(channel_types)
        channel_types = {channel_types};
    end
end

kept_ichans = ismember({channel_def.Channel.Type}, channel_types);
fchannel_def = channel_def;
fchannel_def.Channel = channel_def.Channel(kept_ichans);
fdata = data(:, kept_ichans);

other_ichans = ~kept_ichans;
channel_def_other = channel_def;
channel_def_other.Channel = channel_def.Channel(other_ichans);
data_other = data(:, other_ichans);
end

function [data_concat, channel_def_concat] = ...
    concatenate_data(data, channel_def, data_supp, channel_def_supp)
    
data_concat = [data data_supp];
channel_def_concat = channel_def;
channel_def_concat.Channel = [channel_def_concat.Channel channel_def_supp.Channel];
end


function [nirs_hb, channel_hb_def] = ...
    Compute(nirs_sig, channel_def, age, normalize_method, do_plp_corr, pvf,do_SuperficalRegression,supf_chann)
%% Apply MBLL to compute [HbO] & [HbR] from given nirs OD data
% Args
%    - nirs_sig: matrix of double, size: nb_samples x nb_channels
%        Measured nirs signal (Optical density). Should not contain
%        negative values
%    - channel_def: struct
%        Defintion of channels as given by brainstorm
%        Used fields: Nirs.Wavelengths, Channel
%        ASSUME: channel coordinates are in meters
%    [- age ]: positive double, default is 25
%        Age of the subject, used for light path length correction
%    [- normalize_method]: str in {'mean','median'}, default: 'mean'
%        Method to compute delta optical density values: how to compute the
%        reference intensity against which to compute variations.
%    [- do_ppl_corr]: bool, default: 1
%        Flag to enable partial light path correction (account for light
%        scattering through head tissues)
% 
% Output: 
%   - nirs_hb: matrix of double, size: nb_samples x (nb_channels/nb_wavelengths)*2
%       HbO and HbO delta concentration time-series in mol.l^-1
%   - channel_hb_def: struct 
%       definition of new Hb-related channels. Relevant content:
%           Nirs.Hb = {'HbO', 'HbR', 'HbT'};
%           Channel(ichan1).Name = 'SXDXHbO'; % pair SXDX, HbO component
%           Channel(ichan1).Loc = [C C C]; % pair localization (imported from
%                                          % input channel_def)
%           Channel(ichan1).Group = 'HbO'; 
%           Channel(ichan2).Name = 'SXDXHbR'; % pair SXDX, HbR component
%           Channel(ichan2).Loc = [C C C]; % pair localization (imported from
%                                          % input channel_def, same as paired 
%                                          % channel ichan1)
%           Channel(ichan3).Group = 'HbR'; 
%           Channel(ichan3).Name = 'SXDXHbT'; % pair SXDX, HbT component
%           Channel(ichan3).Loc = [C C C]; % pair localization (imported from
%                                          % input channel_def, same as paired 
%                                          % channel ichan1)
%           Channel(ichan3).Group = 'HbT'; 
% TODO: check negative values
if nargin < 3
    age = 25;
end

if nargin < 4
    normalize_method = 'mean';
end

if nargin < 5
   do_plp_corr = 1; 
end

if nargin < 6
    pvf = 50;
end

[nirs_psig, pair_names, pair_loc, pair_indexes] = group_paired_channels(nirs_sig, channel_def);
pair_distances = cpt_distances(channel_def.Channel, pair_indexes) .* 100; %convert to cm

nb_pairs = length(pair_names);
nb_samples = size(nirs_sig, 1);
nirs_hb_p = zeros(nb_pairs, 3, nb_samples);

if do_SuperficalRegression
    [~,idx_supf_chann] = find(ismember(pair_names,supf_chann{1}'));
    delta_od_supf = zeros(2,nb_samples);
    for i_supf_chann = 1:length(idx_supf_chann)
       delta_od_supf = normalize_nirs(squeeze(nirs_psig(idx_supf_chann(i_supf_chann), :, :)), ...
                              normalize_method)+delta_od_supf; 
    end
    delta_od_supf = delta_od_supf./length(idx_supf_chann);
end
for ipair=1:size(nirs_psig, 1)
    hb_extinctions = nst_get_hb_extinctions(channel_def.Nirs.Wavelengths); % cm^-1.l.mol^-1

    delta_od = normalize_nirs(squeeze(nirs_psig(ipair, :, :)), ...
                              normalize_method);
    
    % SSR regression  Gregg 2010 Frontiers
    if do_SuperficalRegression
        a_data = dot(delta_od',delta_od_supf')./dot(delta_od_supf',delta_od_supf');
        delta_od = delta_od - repmat(a_data,[size(delta_od_supf,2),1])'.*delta_od_supf;
    end
    
    if do_plp_corr
        %TODO: ppf can be computed only once before the loop over pairs
        delta_od_ppf_fixed = fix_ppf(delta_od, channel_def.Nirs.Wavelengths, age, pvf);
    else
        delta_od_ppf_fixed = delta_od;
    end
    
    %mol.l^-1:
	nirs_hb_p(ipair, 1:2, :) = 1 ./ pair_distances(ipair) .* ... % cm^-1
                               pinv(hb_extinctions) * ... mol.cm.l^-1
                               delta_od_ppf_fixed; 
    nirs_hb_p(ipair, 3, :) = sum(nirs_hb_p(ipair, 1:2, :), 2);
end
[nirs_hb, channel_hb_def] = pack_hb_channels(nirs_hb_p, pair_names, pair_loc, channel_def);
end

function [nirs_hb, channel_hb_def] = pack_hb_channels(nirs, pair_names, pair_loc, channel_def_orig)
%% Reshape given nirs data to have size: time x nb_channels
%% Build new Hb-specific channel listing
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
hb_names = {'HbO', 'HbR', 'HbT'};
for ipair=1:size(nirs, 1)
    for ihb=1:length(hb_names)
        nirs_hb(:, ichan) = squeeze(nirs(ipair, ihb, :));
        Channel(ichan).Name = [pair_names{ipair} hb_names{ihb}];
        Channel(ichan).Loc = squeeze(pair_loc(ipair, :, :));
        Channel(ichan).Group = hb_names{ihb};
        Channel(ichan).Comment = [];
        Channel(ichan).Orient = [];
        Channel(ichan).Weight = 1;
        Channel(ichan).Type = 'NIRS';
        ichan = ichan + 1;
    end
end

channel_hb_def.Nirs = rmfield(channel_hb_def.Nirs, 'Wavelengths');

channel_hb_def.Nirs.Hb = hb_names;
channel_hb_def.Channel = Channel;
channel_hb_def.Comment = ['NIRS-BRS sensors (' num2str(length(Channel)) ')'];
end

function delta_od_fixed = fix_ppf(delta_od, wavelengths, age, pvf)
%% Fix given optival density measurements to correct for the differential light path 
%% length: account for light scattering within head tissues
%
% Args:
%     - delta_od: matrix of double, size: nb_wavelengths x time
%         The NIRS OD measurements
%     - wavelengths: array of double, size: nb_wavelengths
%         Wavelengths of the NIRS OD measurements
%    [- age]: double, default is 25
%         The subject's age
%
% Output: matrix of double, size: nb_wavelengths x time
%     Corrected NIRS OD measurements


% Duncan et al 1996:
% dpf = y_0 + a1 * age^a2
dpf_ref_data = [ ...
    [690, 5.38, 0.049, 0.877]; ... % WL, y_0, a1, a2 
    [744, 5.11, 0.106, 0.723]; ... % WL, y_0, a1, a2 
    [807, 4.99, 0.067, 0.814]; ... % WL, y_0, a1, a2 
    [832, 4.67, 0.062, 0.819]; ... % WL, y_0, a1, a2 
    ];

y0 = interp1(dpf_ref_data(:,1), dpf_ref_data(:,2), wavelengths, ...
             'linear', 'extrap');
a1 = interp1(dpf_ref_data(:,1), dpf_ref_data(:,3), wavelengths, ...
             'linear', 'extrap');
a2 = interp1(dpf_ref_data(:,1), dpf_ref_data(:,4), wavelengths, ...
             'linear', 'extrap');
dpf = y0 + a1 .* age.^a2;

% ages = 10:50;
% for ia=1:length(ages)
%     dpfs(:, ia) = y0 + a1 .* ages(ia).^a2;
% end
% plot(ages, dpfs(1, :), 'r'); hold on;
% plot(ages, dpfs(2, :), 'b');

ppf = dpf / pvf;

nb_samples = size(delta_od, 2);
delta_od_fixed = delta_od ./ repmat(ppf', 1, nb_samples);
end

function delta_od = normalize_nirs(nirs_sig, method)
%% Normalize given nirs signal
% Args:
%    - nirs_sig: matrix of double, size:  nb_wavelengths x nb_samples
%        NIRS signal to normalize
%   [- method]: str, choices are: 'mean' and 'median', default is 'mean'
%        Normalization method.
%        * 'mean': divide given nirs signal by its mean, for each wavelength
%        * 'median': divide given nirs signal by its median, for each wavelength
% 
% Output: matrix of double, size: nb_channels x nb_wavelengths
%    Normalized NIRS signal.

switch method
    case 'mean'
        od_ref = mean(nirs_sig, 2);
    case 'median'
        od_ref = median(nirs_sig, 2);
end

nb_samples = size(nirs_sig, 2);

delta_od = -log10( nirs_sig ./ repmat(od_ref, 1, nb_samples) );
end

function distances = cpt_distances(channels, pair_indexes)
% Distance unit is the one of brainstorm (from Channel Loc fields)
% -> meter
distances = zeros(size(pair_indexes, 1), 1);
for ipair=1:size(pair_indexes, 1)
    distances(ipair) = euc_dist(channels(pair_indexes(ipair, 1)).Loc(:,1), ...
                                channels(pair_indexes(ipair, 1)).Loc(:,2));
end
end

function d = euc_dist(p1, p2)
    d = sqrt(sum((p1 - p2).^2));
end


function [paired_nirs, pair_names, pair_loc, pair_indexes] = ...
    group_paired_channels(nirs, channel_def)
%% Reshape the given nirs signal to group paired channels, explode 
%% channel data according to pairs
% Args
%    - nirs_sig: matrix of double, size: nb_samples x nb_channels
%        Nirs signal time-series.
%    - channel_def: struct
%        Defintion of channels as given by brainstorm
%        Used fields: Nirs.Wavelengths, Channel
%
% ASSUME: data contain only wavelength-related channels (no AUX etc.)
%
% TOCHECK WARNING: uses containers.Map which is available with matlab > v2008
%
%  Outputs: 
%     - paired_nirs: array of double, size: nb_pairs x nb_wavelengths xnb_samples
%         Nirs signals, regrouped by pair
%     - pair_names: cell array of str, size: nb_pairs
%         Pair names, format: SXDX
%     - pair_loc: array of double, size: nb_pairs x 3 x 2
%         Pair localization (coordinates of source and detector)
%     - pair_indexes: matrix of double, size: nb_pairs x nb_wavelengths
%         Input channel indexes grouped by pairs
%

nb_wavelengths = length(channel_def.Nirs.Wavelengths);
nb_samples = size(nirs, 1);

pair_to_chans = containers.Map();
for ichan=1:length(channel_def.Channel)
    chan_name = channel_def.Channel(ichan).Name;
    iwl = strfind(chan_name, 'WL'); 
    pair_name = chan_name(1:iwl-1);
    %TODO: keep only channel that are wavelength-related
    wl = str2double(chan_name(iwl+2:end));
    if pair_to_chans.isKey(pair_name)
        wla = pair_to_chans(pair_name);
    else
        wla = zeros(1, nb_wavelengths);
    end
    wla(channel_def.Nirs.Wavelengths==wl) = ichan;
    pair_to_chans(pair_name) = wla;
end
nb_pairs = length(channel_def.Channel) / nb_wavelengths;
pair_names = pair_to_chans.keys;
pair_indexes = zeros(nb_pairs, nb_wavelengths);
pair_loc = zeros(nb_pairs, 3, 2);
paired_nirs = zeros(nb_pairs, nb_wavelengths, nb_samples);
for ipair=1:nb_pairs
    p_indexes = pair_to_chans(pair_names{ipair});
    pair_indexes(ipair, :) = p_indexes;
    paired_nirs(ipair, :, :) = nirs(:, p_indexes)';
    pair_loc(ipair, : , :) = channel_def.Channel(pair_indexes(ipair, 1)).Loc;
end

end

