function varargout = process_nst_sci( varargin )
% process_nst_sci: compute the Scalp Coupling Index
%
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
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
% Authors: Thomas Vincent, 2018
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Scalp Coupling Index';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1201;
sProcess.isSeparator = 0;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Scalp-coupling-index';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','raw'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;


sProcess.options.option_low_cutoff.Comment = 'Bandpass Low cut-off';
sProcess.options.option_low_cutoff.Type    = 'value';
sProcess.options.option_low_cutoff.Value   = {0.5, 'Hz', 4};

sProcess.options.option_high_cutoff.Comment = 'Bandpass High cut-off';
sProcess.options.option_high_cutoff.Type    = 'value';
sProcess.options.option_high_cutoff.Value   = {2.5, 'Hz', 4};

end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInputs)
OutputFiles = {};
for iInput=1:length(sInputs)
    % Load recordings
    if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iInput).FileName);
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
    end
    

    ChannelMat = in_bst_channel(sInputs(iInput).ChannelFile);
    [nirs_ichans, tmp] = channel_find(ChannelMat.Channel, 'NIRS');
    
    signals = sDataIn.F(nirs_ichans,:);
    
    low_cutoff = sProcess.options.option_low_cutoff.Value{1};
    high_cutoff = sProcess.options.option_high_cutoff.Value{1};
    
    fs = 1 / diff(sDataIn.Time(1:2));

    [idx, sci, xpower, xpower_f] = compute(signals,ChannelMat.Channel(nirs_ichans), fs,low_cutoff,high_cutoff );

    % Save time-series data
    data_out = zeros(size(sDataIn.F, 1), length(idx));
    data_out(nirs_ichans,:) = 100*sci;
    sDataOut = db_template('data');
    sDataOut.F            = data_out;
    sDataOut.Comment      = 'SCI';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time(idx);
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.DisplayUnits = '%';
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    OutputFiles{end+1} = OutputFile;

    
    % Save time-series data
    data_out = zeros(size(sDataIn.F, 1), length(idx));
    data_out(nirs_ichans,:) = 100*xpower;
    sDataOut = db_template('data');
    sDataOut.F            = data_out;
    sDataOut.Comment      = 'xpower';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time(idx);
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.DisplayUnits = '%';
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);
    OutputFiles{end+1} = OutputFile;

    % Save time-series data
    data_out = zeros(size(sDataIn.F, 1), length(idx));
    data_out(nirs_ichans,:) = xpower_f;
    sDataOut = db_template('data');
    sDataOut.F            = data_out;
    sDataOut.Comment      = 'Cardiac frequency (xpower)';
    sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
    sDataOut.Time         = sDataIn.Time(idx);
    sDataOut.DataType     = 'recordings';
    sDataOut.nAvg         = 1;
    sDataOut.DisplayUnits = 'Hz';
    
    % Generate a new file name in the same folder
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
    sDataOut.FileName = file_short(OutputFile);
    bst_save(OutputFile, sDataOut, 'v7');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFile, sDataOut);

    OutputFiles{end+1} = OutputFile;
end
end

function [idx_win_start, sci, xpower, xpower_f] = compute(signals, Channel, fs, low_cutoff, high_cutoff, wlen )
%By: Guilherme A. Zimeo Morais 
% https://github.com/GuilhermeZimeo/fNIRS-analysis
%Contact: gazmorais@gmail.com

        % Bandpass filter
     if nargin < 4
         low_cutoff = 0.5;
     end

     if nargin < 5
         high_cutoff = 2.5;
     end

     if nargin < 6
        wlen  = 10; %10s sliding windows 
     end    
     
    order = 3;
    wlen = round(wlen*fs); % convert s to sample
    
    if bst_get('UseSigProcToolbox') 
        [b,a] = butter(order,[low_cutoff*2/fs high_cutoff*2/fs]);
        signals = filtfilt(b, a, signals'); % ntime x nchan
    else
        [b,a] = oc_butter(order,[low_cutoff*2/fs high_cutoff*2/fs]);
        signals = oc_filtfilt(b, a, signals'); % ntime x nchan
    end

    n_sample  = size(signals,1);
    n_channel = size(signals,2);
    
    % Compute 0-lag cross-correlation
    pair_indexes = nst_get_pair_indexes_from_names({Channel(:).Name});

    % windows starts
    idx_win_start = 1:wlen:(n_sample-wlen);

    sci         = zeros(n_channel,length(idx_win_start));
    xpower      = zeros(n_channel,length(idx_win_start));
    xpower_f    = zeros(n_channel,length(idx_win_start));

    for ipair=1:size(pair_indexes, 1)
        chan_indexes = pair_indexes(ipair, 1:2); %Consider only 2 first wavelengths


        %initialize counter for sliding window iteration
        i=1;
        for idx_win=1:length(idx_win_start)

            w = idx_win_start(idx_win);
            
            %normalize data of each wavelength to its standard deviation
            wl1_norm = signals(w:w+wlen,chan_indexes(1))./std(signals(w:w+wlen,chan_indexes(1)),1);
            wl2_norm = signals(w:w+wlen,chan_indexes(2))./std(signals(w:w+wlen,chan_indexes(2)),1);

            
            %calculate scalp coupling index
            [xcorr_wl, lags] = xcorr(wl1_norm,wl2_norm,'coeff');
            sci(chan_indexes,i) = xcorr_wl(lags == 0);

            %Compute Power Spectrum of normalized cross correlation
            L   = length(lags);
            y   = fft(xcorr_wl);
            P2 = abs(y / L );
            P1 = P2(1:round(length(lags)/2));
            P1(2:end-1) = 2*P1(2:end-1);

            freq_fft = fs/L*(0:(L/2));

            %calculate xpower and frequency
            [max_power, idx] = max(P1);

            xpower(chan_indexes,i) = max_power;   
            xpower_f(chan_indexes,i) = freq_fft(idx);   

            i= i+1;
        end

         %final SCI and xpower values (median from sliding window)
        %sci(chan_indexes) = median(SCI_w);
        %xpower(chan_indexes) = median(xpower_w);
    end

end


function sig = normalize(sig)
nb_samples = size(sig, 2);

% sig = 2 * (sig - repmat(min(sig,[],2), 1, nb_samples)) ./ ...
%       repmat(max(sig,[],2)-min(sig,[],2), 1, nb_samples) - 1;
%sig = sig - repmat(mean(sig,2), 1, nb_samples);
sig = sig ./ repmat(std(sig, 0, 2), 1, nb_samples);
end

