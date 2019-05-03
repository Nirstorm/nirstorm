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
% Authors: Thomas Vincent, 2017
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
    sProcess.Comment     = 'Import from CSV table file';
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
        {{'.csv','.xls', '.xpd'}, 'CSV table (*.csv, *.xls, *.xpd)', 'CSV'}, ... % Get all the available file formats
        'EventsIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
    % Option: Event file
    sProcess.options.evtfile.Comment = 'Event file: ';
    sProcess.options.evtfile.Type    = 'filename';
    sProcess.options.evtfile.Value   = SelectOptions;
    
    sProcess.options.label_cols.Comment = '<HTML> <B> Column specifications </B>';
    sProcess.options.label_cols.Type    = 'label';
    
    sProcess.options.trial_label_column.Comment = 'Column name for trial label: ';
    sProcess.options.trial_label_column.Type = 'text';
    sProcess.options.trial_label_column.Value = '';
    
    % Option onset column
    sProcess.options.trial_start_column.Comment = 'Column name for trial start: ';
    sProcess.options.trial_start_column.Type = 'text';
    sProcess.options.trial_start_column.Value = '';
    
    sProcess.options.label_span_help.Comment = 'Trial span:';
    sProcess.options.label_span_help.Type    = 'label';
    [ts, sProcess.options.span_type] = get_trial_spans_opt();
    
    sProcess.options.trial_end_column.Comment = 'Column name for trial end or duration: ';
    sProcess.options.trial_end_column.Type = 'text';
    sProcess.options.trial_end_column.Value = '';
  
    [tu, sProcess.options.time_unit] = get_time_units_opt();
    
    sProcess.options.label_timing.Comment = '<HTML> <B> Time origin </B>';
    sProcess.options.label_timing.Type    = 'label';
  
    [to, sProcess.options.time_origin_type] = get_time_origin_opt();
    
    sProcess.options.time_origin_value_sec.Comment = 'Time origin -- specific value: ';
    sProcess.options.time_origin_value_sec.Type = 'value';
    sProcess.options.time_origin_value_sec.Value = {0.0, 'sec.', 4};
    
    sProcess.options.time_origin_offset_sec.Comment = 'Time origin -- offset: ';
    sProcess.options.time_origin_offset_sec.Type = 'value';
    sProcess.options.time_origin_offset_sec.Value = {0.0, 'sec.', 4};
    
    sProcess.options.label_extra.Comment = '<HTML> <B> Extra options </B>';
    sProcess.options.label_extra.Type    = 'label';
    
    sProcess.options.label2.Comment = 'Entry filters (col_name1=valueA, col_name2<=valueB ...):';
    sProcess.options.label2.Type    = 'label';
    sProcess.options.entry_filters.Comment = '';
    sProcess.options.entry_filters.Type = 'text';
    sProcess.options.entry_filters.Value = '';
    
    sProcess.options.delimiter.Comment = 'Delimiter: ';
    sProcess.options.delimiter.Type = 'text';
    sProcess.options.delimiter.Value = ',';
    
    sProcess.options.max_events.Comment = 'Maximum number of events: ';
    sProcess.options.max_events.Type = 'value';
    sProcess.options.max_events.Value = {1000, '', 0};
    
    sProcess.options.confirm_importation.Comment = 'Confirm importation';
    sProcess.options.confirm_importation.Type = 'checkbox';
    sProcess.options.confirm_importation.Value = 1; 
end



function [time_units, process_option] = get_time_units_opt()
time_units.SECOND = 1;
time_units.MILLISECOND = 2;
time_units.MICROSECOND = 3;
time_units.NANOSECOND = 4;

process_option.Comment = {'sec.', 'msec.', 'musec.', 'nsec.', 'Time unit:'};
process_option.Type    =  'radio_line';
process_option.Value   = time_units.SECOND;
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
    time_units = process_nst_import_csv_events('get_time_units_opt');
    time_origin_types = process_nst_import_csv_events('get_time_origin_opt');
    
    %% Process inputs & do checks
    event_file  = sProcess.options.evtfile.Value{1};
    if ~exist(event_file, 'file')
        bst_error(['Event file not found: ' event_file]);
        return;
    end
    
    event_table = readtable(event_file, 'FileType', 'text', ...
                            'Delimiter', sProcess.options.delimiter.Value, ...
                            'ReadVariableNames', true);
    if size(event_table, 2) < 2
        bst_error(['<HTML> Too few columns (' num2str(size(event_table, 2)) ...
                   ') after loading "' event_file '". <BR> At least 2 columns needed.' ...
                   'Maybe an issue with delimiter or bad CSV format']);
        return;
    end
                        
    if ~check_colname_option(sProcess.options.trial_start_column, event_table)
        return;
    end
    start_colname = sProcess.options.trial_start_column.Value;
    
    if ~check_colname_option(sProcess.options.trial_label_column, event_table)
        return;
    end
    label_colname = sProcess.options.trial_label_column.Value;
    

    
    if sProcess.options.span_type.Value ~= trial_span_types.START_ONLY
        if ~check_colname_option(sProcess.options.trial_end_column, event_table)
            return;
        end
        end_colname = sProcess.options.trial_end_column.Value;
    end
    
    entry_filters = unformat_filters(sProcess.options.entry_filters.Value, event_table);
    if isempty(entry_filters)
        return;
    end
    
%     % Load recordings
%     if strcmp(sInputs.FileType, 'data')     % Imported data structure
%         sDataIn = in_bst_data(sInputs(1).FileName);
%     elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
%         sDataIn = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
%     end
    
    %% Extract events
    selection = apply_filters(entry_filters, event_table);

    if sum(selection) > sProcess.options.max_events.Value{1}
        bst_error('Too many events');
        return
    elseif sum(selection) == 0
        bst_error('No event to import');
        return
    end
    
    switch(sProcess.options.time_unit.Value)
        case time_units.SECOND
            tfactor = 1.;
        case time_units.MILLISECOND
            tfactor = 1e-3;
        case time_units.MICROSECOND
            tfactor = 1e-6;
        case time_units.NANOSECOND
            tfactor = 1e-9;    
    end
    
    event_names = event_table.(label_colname)(selection);
    event_onsets = event_table.(start_colname)(selection);
    if iscell(event_onsets)
        event_onsets = cellfun(@(x) str2double(x), event_onsets) .* tfactor;
    elseif isnumeric(event_onsets)
        event_onsets = event_onsets .* tfactor;
    else
        bst_error(['Column "' start_colname '" could not be converted to numerical values']);
        return
    end
    if size(event_onsets, 1) > 1
        event_onsets = event_onsets';
    end
    if sProcess.options.span_type.Value ~= trial_span_types.START_ONLY
        end_events = event_table.(end_colname)(selection);
        if iscell(end_events)
            end_events = cellfun(@(x) str2double(x), end_events) .* tfactor;
        elseif isnumeric(end_events)
            end_events = end_events .* tfactor;
        else
            bst_error(['Column "' end_colname '" could not be converted to numerical values']);
            return
        end
        if size(end_events, 1) > 1
            end_events = end_events';
        end
        if sProcess.options.span_type.Value == trial_span_types.START_DURATION % end colname is duration
            end_events  = event_onsets + end_events;
        end
    end
    %% Manage time origin
    if sProcess.options.time_origin_type.Value == time_origin_types.OFFSET
        offset_sec = sProcess.options.time_origin_offset_sec.Value{1};
    else % time_origin_types.VALUE
        sorted_onsets = sort(event_onsets);
        offset_sec = sProcess.options.time_origin_value_sec.Value{1} - sorted_onsets(1);
    end
    event_onsets = event_onsets + offset_sec;
    if sProcess.options.span_type.Value ~= trial_span_types.START_ONLY
        end_events  = end_events + offset_sec;
    end

    %% Regroup by trial type
    condition_names = unique(event_names);
    newEvents = repmat(db_template('event'), [1, length(condition_names)]);
    for icond=1:length(condition_names)
        cond_selection = strcmp(event_names, condition_names{icond});
        newEvents(icond).label = condition_names{icond};
        newEvents(icond).times = event_onsets(cond_selection);
        newEvents(icond).epochs = ones(1, length(newEvents(icond).times));
        if sProcess.options.span_type.Value ~= trial_span_types.START_ONLY
            newEvents(icond).times(2,:) = end_events(cond_selection);
        end
    end
        
    if sProcess.options.confirm_importation.Value
        message = {'The following events will be imported:'};
        for ievt=1:length(newEvents)
            if sProcess.options.span_type.Value == trial_span_types.START_ONLY
                duration_info = ' (single)';
            else
                duration_info = sprintf(', avg duration=%1.3f sec.', mean(newEvents(ievt).times(2, :) - newEvents(ievt).times(1, :)));
            end
            message{end+1} = sprintf(' - %s : %d trials. 1st trial at %1.3f sec.%s', ...
                                     newEvents(ievt).label, length(newEvents(ievt).epochs), ...
                                     newEvents(ievt).times(1,1), duration_info);

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

function OutputFile = import_events(sProcess, sInput, newEvents)

OutputFile = [];

% Load the raw file descriptor
isRaw = strcmpi(sInput.FileType, 'raw');
if isRaw
    DataMat = in_bst_data(sInput.FileName, 'F');
    sFile = DataMat.F;
else
    sFile = in_fopen(sInput.FileName, 'BST-DATA');
end

%% ===== MERGE EVENTS LISTS =====
% Add each new event
for iNew = 1:length(newEvents)
    % Look for an existing event
    if ~isempty(sFile.events)
        iEvt = find(strcmpi(newEvents(iNew).label, {sFile.events.label}));
    else
        iEvt = [];
    end
    % Make sure that the sample indices are round values
    if ~isempty(newEvents(iNew).samples)
        newEvents(iNew).times   = newEvents(iNew).samples ./ sFile.prop.sfreq;
    else
        newEvents(iNew).times   = newEvents(iNew).samples ./ sFile.prop.sfreq;
    end
    % If event does not exist yet: add it at the end of the list
    if isempty(iEvt)
        if isempty(sFile.events)
            iEvt = 1;
            sFile.events = newEvents(iNew);
        else
            iEvt = length(sFile.events) + 1;
            sFile.events(iEvt) = newEvents(iNew);
        end
        % Event exists: merge occurrences
    else
        % Merge events occurrences
        %TODO, tell FT that if events are extended and were previously
        % single then bug -> FIX IT
        sFile.events(iEvt).times      = [sFile.events(iEvt).times, newEvents(iNew).times];
        sFile.events(iEvt).epochs     = [sFile.events(iEvt).epochs, newEvents(iNew).epochs];
        sFile.events(iEvt).reactTimes = [sFile.events(iEvt).reactTimes, newEvents(iNew).reactTimes];
        % Sort by sample indices
        if (size(sFile.events(iEvt).samples, 2) > 1)
            [tmp__, iSort] = unique(sFile.events(iEvt).samples(1,:));
            sFile.events(iEvt).times   = sFile.events(iEvt).times(:,iSort);
            sFile.events(iEvt).epochs  = sFile.events(iEvt).epochs(iSort);
            if ~isempty(sFile.events(iEvt).reactTimes)
                sFile.events(iEvt).reactTimes = sFile.events(iEvt).reactTimes(iSort);
            end
        end
    end
    % Add color if does not exist yet
    if isempty(sFile.events(iEvt).color)
        % Get the default color for this new event
        sFile.events(iEvt).color = panel_record('GetNewEventColor', iEvt, sFile.events);
    end
end


%% Set and save the events
if ~isempty(newEvents)
    % Report changes in .mat structure
    if isRaw
        DataMat.F = sFile;
    else
        DataMat.Events = sFile.events;
    end
    % Save file definition
    bst_save(file_fullpath(sInput.FileName), DataMat, 'v6', 1);
    % Report number of detected events
    bst_report('Info', sProcess, sInput, ...
        sprintf('Added to file: %d events in %d different categories', ...
        size([newEvents.epochs],2), length(newEvents)));
    OutputFile = sInput.FileName;
else
    bst_report('Error', sProcess, sInput, 'No events read from file.');
end
end


function selection = apply_filters(entry_filters, event_table)
selection = ones(size(event_table, 1),1) == 1;
if ~ischar(entry_filters) && ~strcmp(entry_filters, 'none')
    for ifilter=1:length(entry_filters)
        filter = entry_filters(ifilter);
        col_to_filter = event_table.(filter.colname);
        
        if strcmp(filter.operator, '=')
            if isnumeric(col_to_filter)
                cur_selection = col_to_filter == str2double(filter.value);
            else
                cur_selection = strcmp(col_to_filter, filter.value);
            end
        elseif strcmp(filter.operator, '~=')
            if isnumeric(col_to_filter)
                cur_selection = col_to_filter ~= str2double(filter.value);
            else
                cur_selection = ~strcmp(col_to_filter, filter.value);
            end
        else
            if ~isnumeric(col_to_filter)
                col_to_filter = cellfun(@(x) str2double(x), col_to_filter);
            end
            filter_value = str2double(filter.value);
            switch(filter.operator)
                case '<'
                    cur_selection = col_to_filter < filter_value;
                case '<='
                    cur_selection = col_to_filter <= filter_value;
                case '>'
                    cur_selection = col_to_filter > filter_value;
                case '>='
                    cur_selection = col_to_filter >= filter_value;
                otherwise
                    bst_error(['Bad entry filter operator: ' filter.operator]);
                    selection = ones(size(event_table, 1),1) == 0;
                    return;
            end
        end
        selection = selection & cur_selection;
    end
    % elseif NIRS start timestamp is defined
    %   -> select trials which occured during NIRS acquisition
end
end

function filters = unformat_filters(sfilters, event_table)

filters = 'none';
if isempty(sfilters)
    return;
end
filter_re = '(([^\s,=<>~]+[~=><]+[^\s,=<>~]+)(?:,\s*)?)+';
match = regexp(sfilters, filter_re, 'match');

if ~iscell(match) || isempty(match) || (~isempty(match) && length(match{1}) ~= length(sfilters))
    bst_error('Malformed entry filter, must be a coma-separated list: "<col_name1><operator><value_A>');
    filters = [];
    return;
end

if ~isempty(sfilters)
    filters = struct([]);
    filter_toks = strsplit(sfilters, ',');
    for ifilter=1:length(filter_toks)
        subtokens = regexp(filter_toks{ifilter}, '([^\s,=<>~]+)([~=><]+)([^\s,=<>~]+)', 'tokens');
        ftoks = subtokens{1};
        filters(ifilter).colname = ftoks{1};
        if ~ismember(ftoks{1}, event_table.Properties.VariableNames)
            bst_error(['Bad entry filter: column "' ftoks{1} '" not found in CSV table.']);
            filters = [];
            return;
        end
        filters(ifilter).operator = ftoks{2};
        if ~ismember(ftoks{2}, {'=', '~=', '<', '<=', '>', '>='})
            bst_error(['Bad entry filter: operator "'  ftoks{2} '".']);
            filters = [];
            return;
        end
        if any(strcmp(filters(ifilter).operator, {'=', '~='}))
            all_values = unique(event_table.(ftoks{1}));
            if isnumeric(all_values)
                value = str2double(ftoks{3}); %TODO: safe conversion
                if isnan(value)
                    bst_error(['Bad entry filter: value "'  ftoks{3} '" is not a number.']);
                    filters = [];
                    return;
                end
            else
                value = ftoks{3};
            end
            if ~ismember(value, all_values)
                bst_error(['Bad entry filter: value "'  ftoks{3} '" not found in column "' ftoks{1} '.']);
                filters = [];
                return;
            end
        end
        filters(ifilter).value = ftoks{3};
    end
end
end


function success = check_colname_option(col_option, event_table)
success = 1;

if isempty(col_option.Value)
    bst_error(['Error. ' col_option.Comment ' not set']);
    success = 0;
    return
end

if ~ismember(col_option.Value, event_table.Properties.VariableNames)
    bst_error(['<HTML>Error. ' col_option.Comment '"' col_option.Value '" is not valid.'...
               ' Available names: <BR>' strjoin(event_table.Properties.VariableNames, ', ')]);
    success = 0;
    return
end
end
