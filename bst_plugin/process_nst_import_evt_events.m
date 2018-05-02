function varargout = process_nst_import_csv_events( varargin )
% PROCESS_EVT_IMPORT: Import events into a raw file

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
% Authors: Edouard Delaire, 2018
%
% Import events from a CSV file with trial as rows and which must at least
% has a column with trial labels and a column with trial onsets. 
% The trial end is optional. If not given, then single events are
% created. If given, it can either be trial durations of trial ending times 
% (the user has to specify this).
%
% The time unit within the CSV file can be specified by the user (sec, msec, musec or nsec).
%
% The time origin of trials can be redefined in two ways:
%    - by adding an offset (in sec.), eg when stimulations and recordings are
%      synchronized on the same time clock and the origin of the stimulations 
%      need to be set relatively to the begining of recordings.
%    - by forcing the origin to a given value (in sec.), eg when stimulations and 
%      recordings were not synchronized but the temporal origin of
%      stimulation has been manually recorded.
%
% Trials, ie rows in the CSV file, can be optionaly filtered by a list of 
% rules. The filter "col_name1=valueA, col_name2<valueB" will select
% only rows where values in col_name1 are equal to valueA and values in
% col_name2 are below valueB.
% Allowed operators: 
%   - "=, ~=": comparison between strings is supported. 
%              If target column can be converted to numerical values then
%              comparison will be numerical (compared value will also be converted).
%   - "<, <=, >, >=": order comparison between numerical values.
%              Target column and value will be converted to numerical
%              values.
% 
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Import from evt table file';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Events';
    sProcess.Index       = 40;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/EventMarkers#Other_menus';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % File selection options
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'open', ...                           % Dialog type: {open,save}
        'Import events...', ...               % Window title
        'ImportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                         % Selection mode: {single,multiple}
        'files', ...                          % Selection mode: {files,dirs,files_and_dirs}
        {{'.evt'},'evt'}, ... % Get all the available file formats
        'EventsIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
    % Option: Event file
    sProcess.options.evtfile.Comment = 'Event file: ';
    sProcess.options.evtfile.Type    = 'filename';
    sProcess.options.evtfile.Value   = SelectOptions;
    
    sProcess.options.label_cols.Comment = '<HTML> <B> Events Names separated by coma </B>';
    sProcess.options.label_cols.Type    = 'text';
    sProcess.options.label_cols.Value    = '';

%     sProcess.options.trial_label_column.Comment = 'Column name for trial label: ';
%     sProcess.options.trial_label_column.Type = 'text';
%     sProcess.options.trial_label_column.Value = '';
    
    
    sProcess.options.label_span_help.Comment = 'Trial span:';
    sProcess.options.label_span_help.Type    = 'label';
    [ts, sProcess.options.span_type] = get_trial_spans_opt();
    
   
    sProcess.options.label_timing.Comment = '<HTML> <B> Time origin </B>';
    sProcess.options.label_timing.Type    = 'label';
  
    [to, sProcess.options.time_origin_type] = get_time_origin_opt();
    
    sProcess.options.time_origin_value_sec.Comment = 'Time origin -- specific value: ';
    sProcess.options.time_origin_value_sec.Type = 'value';
    sProcess.options.time_origin_value_sec.Value = {0.0, 'sec.', 4};
    
    sProcess.options.time_origin_offset_sec.Comment = 'Time origin -- offset: ';
    sProcess.options.time_origin_offset_sec.Type = 'value';
    sProcess.options.time_origin_offset_sec.Value = {0.0, 'sec.', 4};
    
    
    sProcess.options.confirm_importation.Comment = 'Confirm importation';
    sProcess.options.confirm_importation.Type = 'checkbox';
    sProcess.options.confirm_importation.Value = 1; 
end

function [trials_spans, process_option] = get_trial_spans_opt()
trials_spans.START_ONLY = 1;
trials_spans.START_END = 2;
trials_spans.START_DURATION = 3;

process_option.Comment = {'start only  <HTML> <I> (single events, column name for trial end ignored)</I>', ...
                          'start & end  <HTML> <I> (extended events) </I>', ...
                          'start & duration <HTML> <I> (extended events) </I>'};
process_option.Type    = 'radio';
process_option.Value   = trials_spans.START_ONLY;
end


function [time_origin_types, process_option] = get_time_origin_opt()
time_origin_types.OFFSET = 1;
time_origin_types.VALUE = 2;

process_option.Comment = {'Offset', 'Specific', 'Time origin definition:'};
process_option.Type    = 'radio_line';
process_option.Value   = time_origin_types.OFFSET;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};
    
    trial_span_types = process_nst_import_csv_events('get_trial_spans_opt');
    time_origin_types = process_nst_import_csv_events('get_time_origin_opt');
    
    %% Process inputs & do checks
    event_file  = sProcess.options.evtfile.Value{1};
    if ~exist(event_file, 'file')
        bst_error(['Event file not found: ' event_file]);
        return;
    end
    
    event_table = load(event_file);
    n=length(event_table);

    % we create a map for each events events_index(evt) will cointain every
    % time marker related to evt 


    events_names=containers.Map('KeyType','char','ValueType','char');
    events_index=containers.Map('KeyType','char','ValueType','any');

    for i = 1:(n-1)

        key=int2str(bi2de( event_table(i,2:end),'right-msb' ));

        if( isKey (events_index,key) ) 
            events_index(key)=[  events_index(key) ;  event_table(i,1)  event_table(i+1,1)  ];
        else
            events_names(key)=key;
            events_index(key)=[ event_table(i,1)  event_table(i+1,1) ];
        end
    end

    disp('Name #events Mean Duration First Event');
    for key=events_index.keys() 

        evt=events_index(char(key));
        disp( [ key  length(evt) mean( evt(:,2) - evt(:,1) )  evt(1,1) ]);
    end

    
    newEvents = repmat(db_template('event'), [1, length(events_index.keys())]);
    icond=1;
    
    for key=events_index.keys() 
        
        newEvents(icond).label = events_names(char(key));
        newEvents(icond).samples = events_index(char(key))';
        newEvents(icond).epochs = ones(1, length(newEvents(icond).samples));
        
        icond=icond+1;
        
    end
    
    if sProcess.options.confirm_importation.Value
    message = {'The following events will be imported:'};
    for ievt=1:length(newEvents)
        if sProcess.options.span_type.Value == trial_span_types.START_ONLY
            duration_info = ' (single)';
        else
            duration_info = sprintf(', avg duration=%1.3f sec.', mean((newEvents(ievt).samples(2, :) - newEvents(ievt).samples(1, :)))/12.5);
        end
        message{end+1} = sprintf(' - %s : %d trials. 1st trial at %1.3f sec.%s', ...
                                 newEvents(ievt).label, length(newEvents(ievt).epochs), ...
                                 newEvents(ievt).samples(1,1), duration_info);

    end
    
    [res, isCancel] = java_dialog('question', strjoin(message, ['' 10]));
    if strcmp(res, 'No') || isCancel==1
        return;
    end
    end
    
    output_file = import_events(sProcess, sInputs(1), newEvents);

     if ~isempty(output_file)
         OutputFiles{end+1} = output_file;
     end
    
    
end
