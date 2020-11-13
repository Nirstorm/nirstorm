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
sProcess.FileTag     = 'mvt corr';
sProcess.Category    = 'File';
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



sProcess.options.option_event_name.Comment = 'Movement event name: ';
sProcess.options.option_event_name.Type    = 'text';
sProcess.options.option_event_name.Value   = '';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
OutputFile = {};

if ~license('test', 'Curve_Fitting_Toolbox')
    bst_error('Curve Fitting Toolbox not available');
    return 
elseif isempty(which('csaps'))
    bst_error(['Curve Fitting Toolbox OK but function csaps not found.<BR>' ...
               'Try refreshing matlab cache using command: rehash toolboxcache']);
    return
end

% Get selected events
event_name =  strtrim(sProcess.options.option_event_name.Value);


% Load recordings
if strcmp(sInputs.FileType, 'data')     % Imported data structure
    sDataIn = in_bst_data(sInputs.FileName);
    events = sDataIn.Events;
elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file
    sDataIn = in_bst(sInputs.FileName, [], 1, 1, 'no');
    sDataRaw = in_bst_data(sInputs.FileName, 'F');
    events = sDataRaw.F.events;
end
    
event = [];
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


% Process only NIRS channels
channels = in_bst_channel(sInputs.ChannelFile);
nirs_ichans = channel_find(channels.Channel, 'NIRS');
data_nirs = sDataIn.F(nirs_ichans, :)';

prev_negs = any(data_nirs <= 0, 1);

data_corr = Compute(data_nirs, sDataIn.Time', event);

new_negs = any(data_corr <= 0, 1) & ~prev_negs;
negative_chan=find(new_negs);
pair_indexes = nst_get_pair_indexes_from_names({channels.Channel(nirs_ichans).Name});

bst_report('Warning', sProcess, sInputs, 'Motion correction introduced negative values. Will be fixed by local offset');
if any(new_negs)
    
    for ineg=1:length(negative_chan)
        ipair=find(any(pair_indexes(:, :) == negative_chan(ineg),2));
        offset = 2*abs(min(min(data_corr(:, pair_indexes(ipair, :)))));
        data_corr(:, pair_indexes(ipair, :)) = data_corr(:, pair_indexes(ipair, :)) + offset;
        
    end
    [isrcs, idets, measures, channel_type] = nst_unformat_channels({channels.Channel(pair_indexes(ipair, 1)).Name});
    
    msg=sprintf('S%dD%d corrected with offset: %.2f',isrcs,idets,offset);
    bst_report('Warning', sProcess, sInputs, msg);

end

data_corr_full = sDataIn.F';
data_corr_full(:, nirs_ichans) = data_corr;
data_corr = data_corr_full;

% Save time-series data
sDataOut = db_template('data');
sDataOut.F            = data_corr';
sDataOut.Comment      = 'Motion-corrected';
sDataOut.ChannelFlag  = sDataIn.ChannelFlag;
sDataOut.Time         = sDataIn.Time;
sDataOut.DataType     = 'recordings';
sDataOut.History      = sDataIn.History;
sDataOut = bst_history('add', sDataOut, 'process', sProcess.Comment);

sDataOut.nAvg         = 1;
if ~isempty(ievt_mvt)
    sDataOut.Events       = events([1:(ievt_mvt-1) (ievt_mvt+1):length(events)]);
else
    sDataOut.Events       = events;
end
sDataOut.DisplayUnits = sDataIn.DisplayUnits;
    
% Generate a new file name in the same folder
sStudy = bst_get('Study', sInputs.iStudy);
OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'data_motion_corr');
sDataOut.FileName = file_short(OutputFile);
bst_save(OutputFile, sDataOut, 'v7');
% Register in database
db_add_data(sInputs.iStudy, OutputFile, sDataOut);
end


%% ===== Compute =====
function [data_corr] = Compute(nirs_sig, t, event, method) %#ok<DEFNU>
if nargin < 4
    method = 'spline';
end
data_corr = nirs_sig; 
if strcmp(method,'spline') && ~isempty(event)  && ~isempty(event.times)
    samples = time_to_sample_idx(event.times, t);
    data_corr = nst_spline_correction(nirs_sig, t, samples);
end
end

function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end