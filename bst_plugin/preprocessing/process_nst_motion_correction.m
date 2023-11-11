function varargout = process_nst_motion_correction( varargin )

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
% Authors: Thomas Vincent (2015-2018)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
%TOCHECK: how do we limit the input file types (only NIRS data)?
sProcess.Comment     = 'Motion correction';
sProcess.FileTag     = @GetFileTag;
sProcess.Category    = 'Filter';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1305; %0: not shown, >0: defines place in the list of processes
sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Movement_correction';
sProcess.isSeparator = 0; % add a horizontal bar after the process in
%                             the list
% Definition of the input accepted by this process
sProcess.InputTypes  = {'data', 'raw'};
sProcess.OutputTypes = {'data', 'data'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
% Definition of the options


sProcess.options.method.Type       = 'radio_linelabel';
sProcess.options.method.Comment    = {'Spline correction', ' Temporal Derivative Distribution Repair','Motion correction algorithm'; 'spline', 'tddr',''};
sProcess.options.method.Controller = struct('spline','spline','tddr','tddr');
sProcess.options.method.Value      = 'spline';


sProcess.options.option_event_name.Comment = 'Movement event name: ';
sProcess.options.option_event_name.Type    = 'text';
sProcess.options.option_event_name.Value   = '';
sProcess.options.option_event_name.Class   = 'spline';

sProcess.options.option_smoothing.Comment = 'Smoothing Parameters';
sProcess.options.option_smoothing.Type    = 'value';
sProcess.options.option_smoothing.Value   = {0.99,'',3};
sProcess.options.option_smoothing.Class   = 'spline';

sProcess.options.citation.Comment = '<b>Source:</b>';
sProcess.options.citation.Type    = 'label';

sProcess.options.citation_spline.Comment   = ['<p>Scholkmann, F., Spichtig, S., Muehlemann, T., & Wolf, M. (2010). <br />' ...
                                              'How to detect and reduce movement artifacts in near-infrared imaging <br />' ...
                                              'using moving standard deviation and spline interpolation. <br />' ...
                                              'Physiological measurement, 31(5) <br />' ...
                                              'https://doi.org/10.1088/0967-3334/31/5/004<p>'];

sProcess.options.citation_spline.Type    = 'label';
sProcess.options.citation_spline.Class   = 'spline';

sProcess.options.citation_tddr.Comment   = ['<p>Fishburn F.A., Ludlum R.S., Vaidya C.J., & Medvedev A.V. (2019). <br />' ...
                                            'Temporal Derivative Distribution Repair (TDDR): A motion correction <br />' ...
                                             'method for fNIRS. NeuroImage, 184, 171-179. <br />' ...
                                             'https://doi.org/10.1016/j.neuroimage.2018.09.025</p>']; 
sProcess.options.citation_tddr.Type    = 'label';
sProcess.options.citation_tddr.Class    = 'tddr';
end

%% ===== FORMAT COMMENT =====
function [Comment, fileTag] = FormatComment(sProcess)
    % Get options

    % Format comment
     if strcmp(sProcess.options.method.Value,'spline')
        Comment = 'Motion Corrected (spline)';
        fileTag = 'motion';
     else
        Comment = 'Motion Corrected (TDDR)';
        fileTag = 'motion';
    end
end

%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
    [Comment, fileTag] = FormatComment(sProcess);
end


%% ===== RUN =====
function sInputs = Run(sProcess, sInputs) %#ok<DEFNU>

if strcmp(sProcess.options.method.Value,'spline')
    if ~license('test', 'Curve_Fitting_Toolbox')
        bst_error('Curve Fitting Toolbox not available');
        return 
    elseif isempty(which('csaps'))
        bst_error(['Curve Fitting Toolbox OK but function csaps not found.<BR>' ...
                   'Try refreshing matlab cache using command: rehash toolboxcache']);
        return
    end
end

% Get selected events
event_name  =  strtrim(sProcess.options.option_event_name.Value);
events      = sInputs.Events;

% Load recordings
% if strcmp(sInputs.FileType, 'data')     % Imported data structure
%     sDataIn = in_bst_data(sInputs.FileName);
%     events = sDataIn.Events;
% elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
%     sDataIn = in_bst(sInputs.FileName, [], 1, 1, 'no');
%     sDataRaw = in_bst_data(sInputs.FileName, 'F');
%     events = sDataRaw.F.events;
% end
   
event = [];
if strcmp(sProcess.options.method.Value,'spline')
    ievt_mvt = [];
    for ievt=1:length(events)
        if strcmp(events(ievt).label, event_name)
            event = events(ievt);
            ievt_mvt = ievt;
            break;
        end
    end
    if isempty(event)
        warning(['Event "' event_name '" does not exist in file.']);
    end
end



% Process only NIRS channels
channels = in_bst_channel(sInputs.ChannelFile);
nirs_ichans = channel_find(channels.Channel, 'NIRS');
data_nirs = sInputs.A(nirs_ichans, :)';

prev_negs = any(data_nirs <= 0, 1);

data_corr = Compute(data_nirs, sInputs.TimeVector', event,sProcess.options.method.Value,sProcess.options.option_smoothing.Value{1});

new_negs = any(data_corr <= 0, 1) & ~prev_negs;
negative_chan=find(new_negs);
pair_indexes = nst_get_pair_indexes_from_names({channels.Channel(nirs_ichans).Name});

if any(new_negs)
    bst_report('Warning', sProcess, sInputs, 'Motion correction introduced negative values. Will be fixed by local offset');

    for ineg=1:length(negative_chan)
        ipair=find(any(pair_indexes(:, :) == negative_chan(ineg),2));
        offset = 2*abs(min(min(data_corr(:, pair_indexes(ipair, :)))));
        data_corr(:, pair_indexes(ipair, :)) = data_corr(:, pair_indexes(ipair, :)) + offset;
        
    end
    [isrcs, idets, measures, channel_type] = nst_unformat_channels({channels.Channel(pair_indexes(ipair, 1)).Name});
    
    msg=sprintf('S%dD%d corrected with offset: %.2f',isrcs,idets,offset);
    bst_report('Warning', sProcess, sInputs, msg);

end

% Export 
sInputs.A(nirs_ichans,:) = data_corr';
sInputs.CommentTag       = FormatComment(sProcess);

end


%% ===== Compute =====
function [data_corr] = Compute(nirs_sig, t, event, method,exta_parameters) %#ok<DEFNU>
if nargin < 4
    method = 'spline';
end

if nargin < 5
    exta_parameters=0.99;
end    

data_corr = nirs_sig; 
if strcmp(method,'spline') && ~isempty(event)  && ~isempty(event.times)
    samples = time_to_sample_idx(event.times, t);
    data_corr = nst_spline_correction(nirs_sig, t, samples',exta_parameters);
elseif strcmp(method,'tddr')
    fs = 1/(t(2)-t(1));
    data_corr = nst_tddr_correction( nirs_sig , fs );
end
end

function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end