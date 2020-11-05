function varargout = process_nst_iir_filter( varargin )

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
    %TOCHECK: how do we limit the input file types (only NIRS data)?
    sProcess.Comment     = 'Band-pass filter';
    sProcess.FileTag     = @GetFileTag; 
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = {'NIRS', 'Pre-process'};
    sProcess.Index       = 1306; %0: not shown, >0: defines place in the list of processes
    sProcess.isSeparator = 0;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSDataProcess';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'}; %TODO: check processing of link to raw data
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'data'}; %TODO: 'raw' -> 'raw' or 'raw' -> 'data'?
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the options
    
    % === Channel types
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'NIRS';
    sProcess.options.sensortypes.InputTypes = {'data', 'raw'};
    
    % ==== Parameters 
    sProcess.options.label1.Comment = '<BR><U><B>Filtering parameters</B></U>:';
    sProcess.options.label1.Type    = 'label';
    
    
    sProcess.options.option_filter_type.Comment = 'Filter type';
    sProcess.options.option_filter_type.Type    = 'combobox';
    sProcess.options.option_filter_type.Value   = {1, {'bandpass', 'lowpass', 'highpass'}};
    sProcess.options.option_filter_type.Hidden  = 1; % Maintain compatibility
    
    sProcess.options.option_keep_mean.Comment = 'Keep mean';
    sProcess.options.option_keep_mean.Type    = 'checkbox';
    sProcess.options.option_keep_mean.Value   = 1;
    
    sProcess.options.option_low_cutoff.Comment = 'Lower cutoff frequency (0=disable):';
    sProcess.options.option_low_cutoff.Type    = 'value';
    sProcess.options.option_low_cutoff.Value   = {0.01, '', 4};
    
    sProcess.options.option_high_cutoff.Comment = 'Upper cutoff frequency (0=disable):';
    sProcess.options.option_high_cutoff.Type    = 'value';
    sProcess.options.option_high_cutoff.Value   = {0.5, '', 4}; 
    
    sProcess.options.order.Comment = 'Filter order:';
    sProcess.options.order.Type    = 'value';
    sProcess.options.order.Value   = {3, '', 0};
    
        % === Display properties
    sProcess.options.display.Comment = {'process_nst_iir_filter(''DisplaySpec'',iProcess,sfreq);', '<BR>', 'View filter response'};
    sProcess.options.display.Type    = 'button';
    sProcess.options.display.Value   = [];
    

end

%% ===== GET OPTIONS =====
function [HighPass, LowPass,FilterType,Order,KeepMean] = GetOptions(sProcess)
    HighPass = sProcess.options.option_low_cutoff.Value{1};
    LowPass  = sProcess.options.option_high_cutoff.Value{1};
    Order= sProcess.options.order.Value{1}; 
    ftypes = sProcess.options.option_filter_type.Value{2};
    filter_type = ftypes{sProcess.options.option_filter_type.Value{1}};
    
    FilterType='bandpass';
    
    if (HighPass == 0 || isequal(filter_type,'lowpass')) 
        HighPass = [];
        FilterType='lowpass';
    end
    
    if (LowPass == 0 || isequal(filter_type,'highpass'))  
        LowPass = [];
        FilterType='highpass';
    end
    
    
    KeepMean = sProcess.options.option_keep_mean.Value;
end


%% ===== FORMAT COMMENT =====
function [Comment, fileTag] = FormatComment(sProcess)
    % Get options
    [HighPass, LowPass,FilterType] = GetOptions(sProcess);
    % Format comment
    
    switch FilterType
        case 'bandpass'
            Comment = ['Band-pass: ' num2str(HighPass) 'Hz-' num2str(LowPass) 'Hz'];
            fileTag = '_IIR-band';
        case 'highpass'
            Comment = ['High-pass: ' num2str(HighPass) 'Hz'];
            fileTag = '_IIR-high';
        case 'lowpass'
            Comment = ['Low-pass: ' num2str(LowPass) 'Hz'];
            fileTag = '_IIR-low';
    end
end

%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
    [Comment, fileTag] = FormatComment(sProcess);
end


%% ===== RUN =====
function sInputs = Run(sProcess, sInputs) %#ok<DEFNU>

    warning('Deprecated process: use Pre-process > Band-pass filter instead');

    % Get option values
    
    [low_cutoff, high_cutoff,filter_type,order,keep_mean] = GetOptions(sProcess);
    fs = 1 / diff(sInputs.TimeVector(1:2)); 
    
    %TOCHECK: take care of bad channels here?
    % channels = in_bst_channel(sInputs.ChannelFile);
    [nirs_filtered, transient] = Compute(sInputs.A', fs, filter_type, low_cutoff, ...
                            high_cutoff, order, keep_mean);
    sInputs.A = nirs_filtered';
    
    % Add events to represent the transients (edge effects)
    if transient > 0 
        % Time windows with filter transients (two extended events)
        trans = [sInputs.TimeVector(1),                      sInputs.TimeVector(end) - transient; ...
                 sInputs.TimeVector(1) + transient, sInputs.TimeVector(end)];
        % Create a new event type
        sInputs.Events = db_template('event');
        sInputs.Events.label    = ['transient_' filter_type];
        sInputs.Events.color    = [.8 0 0];
        sInputs.Events.epochs   = [1 1];
        sInputs.Events.times    = round(trans .* fs) ./ fs;
        sInputs.Events.channels = cell(1, size(sInputs.Events.times, 2));
        sInputs.Events.notes    = cell(1, size(sInputs.Events.times, 2));
    end
    
    % File comment
    switch filter_type
        case 'bandpass'
            filterComment = ['band(' num2str(low_cutoff) '-' num2str(high_cutoff) 'Hz)'];
        case  'highpass'   
            filterComment = ['high(' num2str(low_cutoff) 'Hz)'];
        case 'lowpass'
            filterComment = ['low(', num2str(high_cutoff) 'Hz)'];
    end
    sInputs.CommentTag = filterComment;
    
    % Do not keep the Std field in the output
    if isfield(sInputs, 'Std') && ~isempty(sInputs.Std)
        sInputs.Std = [];
    end
    
end


%% ===== Compute =====
function [nirs_filtered, transient] = Compute(nirs_sig, fs, filter_type, low_cutoff, ...
                                  high_cutoff, order, keep_mean)
%% Apply Zero-phase Infinite Impulse filtering (Butterworth) on NIRS data.
%
% Args
%    - nirs_sig: matrix of double, size: time x nb_channels 
%        nirs signals to be filtered
%    - filter_type: str], choice in {'bandpass'} default: 'bandpass'
%        type of filter to apply
%    - low_cutoff: double
%        Lower frequency boundary for the filter
%    - high_cutoff: double
%        Higher frequency boundary for the filter
%    - keep_mean: int
%        If 0, then substract the signal mean (remove flat trend)
%
% Output:
%    - nirs_detrend: matrix of double, size: time x nb_channels 
%        detrended nirs signals

    data = double(nirs_sig);


    % Remove the mean before filtering    
    if keep_mean
        nSamples = size(data,1);
        meanMatrix = repmat(mean(data),nSamples,1);

        data= data-meanMatrix;
    end 
    
    [b,a] = ComputeFilter(fs, filter_type, low_cutoff,high_cutoff, order);
    transient = ComputeTransient(b,a,fs);
    
    nirs_filtered = filtfilt(b,a,data);

    % Add the mean back
    if keep_mean
        nirs_filtered=nirs_filtered+meanMatrix;
    end

end


function [b,a] = ComputeFilter(fs, filter_type, low_cutoff,high_cutoff, order)

    switch filter_type
        % HP filtering
        case 'highpass'
            [b,a]= butter(order,low_cutoff*2/fs, 'high');
        % LP filtering
        case 'lowpass'
            [b,a]= butter(order,high_cutoff*2/fs, 'low');
        % BP filtering
        case 'bandpass'
            [b,a]= butter(order,[low_cutoff*2/fs high_cutoff*2/fs]);
        case 'bandstop'
            [b,a]= butter(order,[low_cutoff*2/fs high_cutoff*2/fs], 'stop');
    end
    %figure(1);
    %zplane(b,a)
end

function transient = ComputeTransient(b,a,sfreq)
% Compute transient --  see bst_bandpass_hfilter.m
% Compute the cumulative energy of the impulse response

    % Compute filter response
    if bst_get('UseSigProcToolbox')
        [Ht,t] = impz(b, a, [], sfreq);
    else
        [Ht,t] = oc_impz( b,  a, [], sfreq);
    end

    E = cumsum(Ht.^2) ;
    E = E ./ max(E) ;
    % Compute the effective transient: Number of samples necessary for having 99% of the impulse response energy
    [tmp, iE99] = min(abs(E - 0.99)) ;
    transient      = iE99 / sfreq ;  
end
%% ===== DISPLAY FILTER SPECS =====
function DisplaySpec(iProcess, sfreq) %#ok<DEFNU>
    % Get current process options
    global GlobalData;
    sProcess = GlobalData.Processes.Current(iProcess);

    % Get options
    [HighPass, LowPass,FilterType,order,KeepMean] = GetOptions(sProcess);
    [b,a] = ComputeFilter(sfreq, FilterType, HighPass,LowPass, order);

    % Compute filter response
    if bst_get('UseSigProcToolbox')
        [Hf,Freqs] = freqz(b, a, 2^14, sfreq);
        [Ht,t] = impz(b, a, [], sfreq);
    else
        [Hf,Freqs] = oc_freqz(b,  a, 2^14, sfreq);
        [Ht,t] = oc_impz( b,  a, [], sfreq);
    end
    

    % Check Stability 
    stable = isstable(b,a);
    
    % Compute transient 
    transient = ComputeTransient(b,a,sfreq);
    
    % Configure 
    if isempty(LowPass) || (LowPass == 0)
        XFreqLim = [0, min(3*HighPass, max(Freqs))] ;
    else
        XFreqLim = [0, min(1.5*LowPass, max(Freqs))] ; 
    end

    % Filter description: Left panel
     strFilter1 = ['<HTML>Zero phase <B>IIR filter</B>' '<BR>'];
     if ~isempty(HighPass) && (HighPass > 0) && ~isempty(LowPass) && (LowPass > 0)
         strFilter1 = [strFilter1 'Band-pass: &nbsp;&nbsp;<B>' num2str(HighPass) '-' num2str(LowPass) ' Hz</B><BR>'];
         %strFilter1 = [strFilter1 'Low transition: &nbsp;&nbsp;<B>' num2str(FiltSpec.fcuts(1)) '-' num2str(FiltSpec.fcuts(2)) ' Hz</B><BR>'];
         %strFilter1 = [strFilter1 'High transition: &nbsp;&nbsp;<B>' num2str(FiltSpec.fcuts(3)) '-' num2str(FiltSpec.fcuts(4)) ' Hz</B><BR>'];
     elseif ~isempty(HighPass) && (HighPass > 0)
         strFilter1 = [strFilter1 'High-pass: &nbsp;&nbsp;<B>' num2str(HighPass) ' Hz</B><BR>'];
         %strFilter1 = [strFilter1 'Transition: &nbsp;&nbsp;<B>' num2str(FiltSpec.fcuts(1)) '-' num2str(FiltSpec.fcuts(2)) ' Hz</B><BR>'];
     elseif ~isempty(LowPass) && (LowPass > 0)
         strFilter1 = [strFilter1 'Low-pass: &nbsp;&nbsp;<B>' num2str(LowPass) ' Hz</B><BR>'];
         %strFilter1 = [strFilter1 'Transition: &nbsp;&nbsp;<B>' num2str(FiltSpec.fcuts(1)) '-' num2str(FiltSpec.fcuts(2)) ' Hz</B><BR>'];
     end
    if KeepMean
         strFilter1 = [strFilter1 'Conserve mean: &nbsp;&nbsp;<B>yes</B><BR>'];
    else
         strFilter1 = [strFilter1 'Conserve mean: &nbsp;&nbsp;<B>no</B><BR>'];
     end

    % Filter description: Right panel
    strFilter2 = '<HTML>';
    strFilter2 = [strFilter2 'Filter type: &nbsp;&nbsp;<B> Butterworth </B><BR>'];
    strFilter2 = [strFilter2 'Filter order: &nbsp;&nbsp;<B>' num2str(order) '</B><BR>'];
    if stable % Maybe add a warning if filter is not stable ? 
        strFilter2 = [strFilter2 'Stable: &nbsp;&nbsp;<B> yes</B><BR>'];
    else
        strFilter2 = [strFilter2 'Stable: &nbsp;&nbsp;<B> no</B><BR>'];
    end    
    strFilter2 = [strFilter2 'Transient (99% energy): &nbsp;&nbsp;<B>' num2str(transient, '%1.3f') ' s</B><BR>'];
    strFilter2 = [strFilter2 'Sampling frequency: &nbsp;&nbsp;<B>', num2str(sfreq), ' Hz</B><BR>'];
    %strFilter2 = [strFilter2 'Frequency resolution: &nbsp;&nbsp;<B>' num2str(dF, '%1.3f') ' Hz</B><BR>'];
    
    hFig = HFilterDisplay(Hf,Freqs,Ht,t,transient,strFilter1,strFilter2,XFreqLim) ; 

end

function hFig = HFilterDisplay(Hf,Freqs,Ht,t,transient,strFilter1,strFilter2,XFreqLim)

% Progress bar
bst_progress('start', 'Filter specifications', 'Updating graphs...');

% Get existing specification figure
hFig = findobj(0, 'Type', 'Figure', 'Tag', 'FilterSpecs');
% If the figure doesn't exist yet: create it
if isempty(hFig)
    hFig = figure(...
        'MenuBar',     'none', ...
        ... 'Toolbar',     'none', ...
        'Toolbar',     'figure', ...
        'NumberTitle', 'off', ...
        'Name',        sprintf('Filter properties'), ...
        'Tag',         'FilterSpecs', ...
        'Units',       'Pixels');
    % Figure already exists: re-use it
else
    clf(hFig);
    figure(hFig);
end

% Plot frequency response
hAxesFreqz = axes('Units', 'pixels', 'Parent', hFig, 'Tag', 'AxesFreqz');
% transfer function equal to the squared magnitude of the original filter
% transfer function. (filtfilt) 


Hf(1)=Hf(2); % better display 

Hf = 40.*log10(abs(Hf));
plot(hAxesFreqz, Freqs, Hf);
% Plot impulse response
hAxesImpz = axes('Units', 'pixels', 'Parent', hFig, 'Tag', 'AxesImpz');
plot(hAxesImpz, t, Ht);

% Add Axes limits
set(hAxesFreqz, 'XLim', XFreqLim);
set(hAxesFreqz, 'YLim', [min(Hf(Freqs<XFreqLim(2))), max(Hf(Freqs<XFreqLim(2)))] + (max(Hf(Freqs<XFreqLim(2)))-min(Hf(Freqs<XFreqLim(2)))) .* [-0.05,0.05]);
YLimImpz = [min(Ht), max(Ht)] + (max(Ht)-min(Ht)) .* [-0.05,0.05];
set(hAxesImpz, 'XLim', [min(t), max(t)], 'YLim', YLimImpz);

% Add grids
set([hAxesFreqz, hAxesImpz], 'XGrid', 'on', 'YGrid', 'on');
% Enable zooming by default
zoom(hFig, 'on');
    
% Add legends
title(hAxesFreqz, 'Frequency response');
xlabel(hAxesFreqz, 'Frequency (Hz)');
ylabel(hAxesFreqz, 'Magnitude (dB)');
title(hAxesImpz, 'Impulse response');
xlabel(hAxesImpz, 'Time (seconds)');
ylabel(hAxesImpz, 'Amplitude');

% Plot vertical lines to indicate effective transients (99% energy)
line(transient.*[1 1], YLimImpz, -0.1.*[1 1], ...
    'LineWidth', 1, ...
    'Color',     [.7 .7 .7], ...
    'Parent',    hAxesImpz);
line(transient.*[-1 -1], YLimImpz, -0.1.*[1 1], ...
    'LineWidth', 1, ...
    'Color',     [.7 .7 .7], ...
    'Parent',    hAxesImpz);
text(transient .* 1.1, YLimImpz(2), '99% energy', ...
    'Color',               [.7 .7 .7], ...
    'FontSize',            bst_get('FigFont'), ...
    'FontUnits',           'points', ...
    'VerticalAlignment',   'top', ...
    'HorizontalAlignment', 'left', ...
    'Parent',              hAxesImpz);

% Display left panel
[jLabel1, hLabel1] = javacomponent(javax.swing.JLabel(strFilter1), [0 0 1 1], hFig);
set(hLabel1, 'Units', 'pixels', 'BackgroundColor', get(hFig, 'Color'), 'Tag', 'Label1');
bgColor = get(hFig, 'Color');
jLabel1.setBackground(java.awt.Color(bgColor(1),bgColor(2),bgColor(3)));
jLabel1.setVerticalAlignment(javax.swing.JLabel.TOP);

% Display right panel
[jLabel2, hLabel2] = javacomponent(javax.swing.JLabel(strFilter2), [0 0 1 1], hFig);
set(hLabel2, 'Units', 'pixels', 'BackgroundColor', get(hFig, 'Color'), 'Tag', 'Label2');
bgColor = get(hFig, 'Color');
jLabel2.setBackground(java.awt.Color(bgColor(1),bgColor(2),bgColor(3)));
jLabel2.setVerticalAlignment(javax.swing.JLabel.TOP);

% Set resize function
set(hFig, bst_get('ResizeFunction'), @ResizeCallback);
% Force calling the resize function at least once
ResizeCallback(hFig);
bst_progress('stop');

% Resize function
    function ResizeCallback(hFig, ev)
        % Get figure position
        figpos = get(hFig, 'Position');
        textH = 110;        % Text Height
        marginL = 70;
        marginR = 30;
        marginT = 30;
        marginB = 50;
        axesH = round((figpos(4) - textH) ./ 2);
        % Position axes
        set(hAxesFreqz, 'Position', max(1, [marginL, textH + marginB + axesH, figpos(3) - marginL - marginR, axesH - marginB - marginT]));
        set(hAxesImpz,  'Position', max(1, [marginL, textH + marginB,         figpos(3) - marginL - marginR, axesH - marginB - marginT]));
        set(hLabel1,    'Position', max(1, [40,                  1,  round((figpos(3)-40)/2),  textH]));
        set(hLabel2,    'Position', max(1, [round(figpos(3)/2),  1,  round(figpos(3)/2),       textH]));
    end
end

