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
sProcess.Comment     = 'Signal Quality Control';
sProcess.Category    = 'File';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1201;
sProcess.isSeparator = 0;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Scalp-coupling-index';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','raw'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;


sProcess.options.window_length.Comment = 'Window length';
sProcess.options.window_length.Type    = 'value';
sProcess.options.window_length.Value   = {10, 's', 0};

sProcess.options.text1.Comment   = '<b>Export the following indicator</b>'; 
sProcess.options.text1.Type    = 'label';


sProcess.options.option_coefficient_variation.Comment = 'Coefficient of variation';
sProcess.options.option_coefficient_variation.Type    = 'checkbox';
sProcess.options.option_coefficient_variation.Value   = 0;
sProcess.options.option_coefficient_variation.Controller='variation';


sProcess.options.option_sci.Comment = 'Scalp Coupling Index';
sProcess.options.option_sci.Type    = 'checkbox';
sProcess.options.option_sci.Value   = 0;
sProcess.options.option_sci.Controller='sci';

sProcess.options.option_low_cutoff.Comment = 'Bandpass Low cut-off';
sProcess.options.option_low_cutoff.Type    = 'value';
sProcess.options.option_low_cutoff.Value   = {0.5, 'Hz', 4};
sProcess.options.option_low_cutoff.Class='sci';

sProcess.options.option_high_cutoff.Comment = 'Bandpass High cut-off';
sProcess.options.option_high_cutoff.Type    = 'value';
sProcess.options.option_high_cutoff.Value   = {2.5, 'Hz', 4};
sProcess.options.option_high_cutoff.Class='sci';


end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
Comment = sProcess.Comment;
end

function OutputFiles = Run(sProcess, sInputs)

    OutputFiles = {};

    % Load recordings
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs.FileName);
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs.FileName, [], 1, 1, 'no');
    end
    
    fs = 1 / diff(sDataIn.Time(1:2));
    window_length = sProcess.options.window_length.Value{1};

    ChannelMat          = in_bst_channel(sInputs.ChannelFile);
    [nirs_ichans, tmp]  = channel_find(ChannelMat.Channel, 'NIRS');
    signals = sDataIn.F(nirs_ichans,:);

    if sProcess.options.option_coefficient_variation.Value

        CV = compute_CV(sDataIn.Time, signals, window_length);

        % Save time-series data
        data_out = zeros(size(sDataIn.F));
        data_out(nirs_ichans,:) = 100 * CV;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'Coefficent of Variation';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = sDataIn.Time;
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = '%';
    
        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputs.iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs.iStudy, OutputFile, sDataOut);
        OutputFiles{end+1} = OutputFile;

    end

    if sProcess.options.option_sci.Value

        low_cutoff = sProcess.options.option_low_cutoff.Value{1};
        high_cutoff = sProcess.options.option_high_cutoff.Value{1};
        

        [idx, sci, xpower, xpower_f] = compute(signals,ChannelMat.Channel(nirs_ichans), fs, low_cutoff, high_cutoff );
        idx = sDataIn.Time;

        % Save time-series data
        data_out = zeros(size(sDataIn.F, 1), length(idx));
        data_out(nirs_ichans,:) = 100*sci;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'SCI';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = sDataIn.Time;
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = '%';
    
        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputs.iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs.iStudy, OutputFile, sDataOut);
        OutputFiles{end+1} = OutputFile;

    
        % Save time-series data
        data_out = zeros(size(sDataIn.F, 1), length(idx));
        data_out(nirs_ichans,:) = 100*xpower;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'xpower';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = sDataIn.Time;
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = '%';
    
        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputs.iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs.iStudy, OutputFile, sDataOut);
        OutputFiles{end+1} = OutputFile;

        % Save time-series data
        data_out = zeros(size(sDataIn.F, 1), length(idx));
        data_out(nirs_ichans,:) = xpower_f;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'Cardiac frequency (xpower)';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = sDataIn.Time;
        sDataOut.DataType     = 'recordings';
        sDataOut.nAvg         = 1;
        sDataOut.DisplayUnits = 'Hz';
    
        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputs.iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_sci');
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs.iStudy, OutputFile, sDataOut);
    
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

    sci         = zeros(n_channel, n_sample);
    xpower      = zeros(n_channel, n_sample);
    xpower_f    = zeros(n_channel, n_sample);

    for ipair=1:size(pair_indexes, 1)
        chan_indexes = pair_indexes(ipair, 1:2); %Consider only 2 first wavelengths


        %initialize counter for sliding window iteration
        for i = 1:n_sample
    
            windows_start =  max( 1, i - wlen/2);
            windows_end   =  min( n_sample, i + wlen/2);
    
            nirs_windows  = signals(windows_start:windows_end, chan_indexes);

            
            %normalize data of each wavelength to its standard deviation
            wl1_norm = nirs_windows(:,1) ./ std(nirs_windows(:,1));
            wl2_norm = nirs_windows(:,2) ./ std(nirs_windows(:,2));

            
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
        end

    end

end


function CV = compute_CV(Time, signals, wlen)
    if nargin < 6
        wlen  = 10; %10s sliding windows 
    end   

    fs = 1 / diff(Time(1:2));
    wlen = round(wlen*fs); % convert s to sample

    n_channel = size(signals,1);
    n_sample = size(signals,2);

    mov_mean    = zeros(n_channel,n_sample); 
    mov_std     = zeros(n_channel,n_sample); 

    for i = 1:n_sample

       windows_start =  max( 1,i - wlen/2);
       windows_end   =  min( n_sample, i + wlen/2);

       nirs_windows  = signals(:, windows_start:windows_end);
       mov_mean(:,i) = mean(nirs_windows, 2 );
       mov_std(:, i) = std(nirs_windows, [], 2);
    end    

    CV = mov_mean ./ mov_std;

end



