function varargout = process_nst_quality_check( varargin )
% process_nst_quality_check: compute several indicator of the signal
% quality
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
% Authors: Edouard Delaire, 2025
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
% Description the process
sProcess.Comment     = 'Signal Quality Control';
sProcess.Category    = 'File';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1201;
sProcess.isSeparator = 0;
sProcess.Description = 'https://neuroimage.usc.edu/brainstorm/Tutorials/NIRSTORM#Signal_quality_check';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','raw'};

sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;

sProcess.options.text0.Comment   = '<b>Options</b>'; 
sProcess.options.text0.Type    = 'label';

sProcess.options.window_length.Comment = 'Window length';
sProcess.options.window_length.Type    = 'value';
sProcess.options.window_length.Value   = {10, 's', 0};

sProcess.options.text1.Comment   = '<b>Export the following indicator</b>'; 
sProcess.options.text1.Type    = 'label';


sProcess.options.option_coefficient_variation.Comment = 'Coefficient of variation';
sProcess.options.option_coefficient_variation.Type    = 'checkbox';
sProcess.options.option_coefficient_variation.Value   = 1;
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

sProcess.options.text2.Comment   = '<b>Export the following figures</b>'; 
sProcess.options.text2.Type    = 'label';

sProcess.options.option_light_intensity.Comment = 'Light intensity';
sProcess.options.option_light_intensity.Type    = 'checkbox';
sProcess.options.option_light_intensity.Value   = 0;
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
        sDataRaw = in_bst_data(sInputs.FileName, 'F');
        sDataIn.Events = sDataRaw.F.events;
    end
    
    window_length = sProcess.options.window_length.Value{1};

    ChannelMat          = in_bst_channel(sInputs.ChannelFile);
    [nirs_ichans, tmp]  = channel_find(ChannelMat.Channel, 'NIRS');
    signals = sDataIn.F(nirs_ichans,:);

    if sProcess.options.option_coefficient_variation.Value

        [Time, CV] = compute_CV(sDataIn.Time, signals, window_length);

        % Save time-series data
        data_out = zeros(size(sDataIn.F));
        data_out(nirs_ichans,:) = CV;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'Coefficent of Variation';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = Time;
        sDataOut.Events       = sDataIn.Events;
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
        

        [Time, sci, xpower, xpower_f] = compute_SCI(sDataIn.Time, signals, window_length, ChannelMat.Channel(nirs_ichans), low_cutoff, high_cutoff );

        % Save time-series data
        data_out = zeros(size(sDataIn.F));
        data_out(nirs_ichans,:) = sci;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'Scalp Coupling Index';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = Time;
        sDataOut.Events       = sDataIn.Events;
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
        data_out = zeros(size(sDataIn.F));
        data_out(nirs_ichans,:) =  xpower;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'Cardiac signal strength';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = Time;
        sDataOut.Events       = sDataIn.Events;
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
        data_out = zeros(size(sDataIn.F));
        data_out(nirs_ichans,:) = xpower_f;
        sDataOut = db_template('data');
        sDataOut.F            = data_out;
        sDataOut.Comment      = 'Cardiac frequency (xpower)';
        sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
        sDataOut.Time         = Time;
        sDataOut.Events       = sDataIn.Events;
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

    if sProcess.options.option_light_intensity.Value 
        plot_intensity(sDataIn.Time, signals, window_length, ChannelMat.Channel(nirs_ichans));
    end
end

function [Time, sci, xpower, xpower_f] = compute_SCI(Time, signals, wlen, Channel, low_cutoff, high_cutoff )
     
    if nargin < 3 
        wlen  = 10; 
     end    

     % Bandpass filter
     if nargin < 5
         low_cutoff = 0.5;
     end

     if nargin < 6
         high_cutoff = 2.5;
     end

    order = 3;

    fs = 1 / diff(Time(1:2));
    wlen = round(wlen*fs);   

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

    sci         = zeros(n_channel, n_sample);
    xpower      = zeros(n_channel, n_sample);
    xpower_f    = zeros(n_channel, n_sample);
    
    L   = 2*wlen + 1;
    freq_fft = fs/L*(0:(L/2));

    bst_progress('start', 'Quality check', 'Computing SCI', 0, size(pair_indexes, 1));

    for ipair=1:size(pair_indexes, 1)

        %Consider only 2 first wavelengths
        chan_indexes = pair_indexes(ipair, 1:2); 

        % Sliding window, with an overlap of 50%
        for i = 1:wlen/2:n_sample
    
            windows_start =  max( 1, i - wlen/2);
            windows_end   =  min( n_sample, i + wlen/2);
    
            nirs_windows  = signals(windows_start:windows_end, chan_indexes);

            
            %normalize data of each wavelength to its standard deviation
            wl1_norm = nirs_windows(:,1) ./ std(nirs_windows(:,1));
            wl2_norm = nirs_windows(:,2) ./ std(nirs_windows(:,2));

            
            %calculate scalp coupling index
            [xcorr_wl, lags]    = xcorr(wl1_norm,wl2_norm,'coeff');
            sci(chan_indexes, windows_start:windows_end) = xcorr_wl(lags == 0);

            y   = fft(xcorr_wl);
            P2 = abs(y ./ L );
            P1 = P2(1:round(length(lags)/2));
            P1(2:end-1) = 2*P1(2:end-1);

            %calculate xpower and frequency
            [max_power, idx] = max(P1);

            xpower(chan_indexes,   windows_start:windows_end) = max_power;   
            xpower_f(chan_indexes, windows_start:windows_end) = freq_fft(idx);  
        end

        bst_progress('inc', 1);
    end

    bst_progress('stop');
end


function [Time, CV] = compute_CV(Time, signals, wlen)
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

    CV =  mov_std ./ mov_mean  ;

end

function [Time, hFig] = plot_intensity(Time, signals, wlen, Channel)
    % Plot the light fall off. Light intensity as function of the
    % source-detector distance
    

    groups = unique({Channel.Group});

    hFig = figure('Color', 'k', 'Position', [200, 200, 560, 420]);
    hold on;

    for iGroup = 1:length(groups)
        idx_chan = strcmp({Channel.Group},groups{iGroup});

        separation_group = process_nst_separations('Compute', Channel(idx_chan)) * 100;
        intensity_group  = mean(signals(idx_chan,:),2);

        semilogy(separation_group, intensity_group,'.');
    end

    intensity  = mean(abs(signals),2);
    M  = ceil(max(log10(abs(intensity(:)))));
    yM = 10^M;

    xM = min([max(separation_group)+0.5,100]);

    axis([0,xM, 1e-6, yM]);
    yscale log
    xlabel('Source-Detector Separation ( mm )','Color','w')
    ylabel('{\Phi_0} ( {\mu}W )','Color','w')
    set(gca,'XColor','w','YColor','w','Xgrid','on','Ygrid','on','Color','k')
    legend(groups,'Color','w')


end

