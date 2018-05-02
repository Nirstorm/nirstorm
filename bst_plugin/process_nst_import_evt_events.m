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
% Import events from a evt file with the following structure 
% sample1  X X X X  
% sample2  Y Y Y Y
%       ... 
% such as the event coded by XXXX last from sample1 to sample2

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

    sProcess.options.confirm_importation.Comment = 'Confirm importation';
    sProcess.options.confirm_importation.Type = 'checkbox';
    sProcess.options.confirm_importation.Value = 1; 
end




%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};
    
    trial_span_types = process_nst_import_csv_events('get_trial_spans_opt');
    
    %% Process inputs & do checks
    event_file  = sProcess.options.evtfile.Value{1};
    raw_events_name= strsplit(sProcess.options.label_cols.Value, ','); 
    number_of_events=length(raw_events_name);
    
    
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
    
    
    k=1;
    
    for i = 1:(n-1)

        key=int2str(bi2de( event_table(i,2:end),'right-msb' ));

        if( isKey (events_index,key) ) 
            events_index(key)= [  events_index(key) ;  event_table(i,1)  event_table(i+1,1)  ];
        else
            if k > number_of_events 
                bst_error('Not enough event name');
                return;
            end
            events_names(key)=char(raw_events_name(k));
            k=k+1;
            events_index(key)= [ event_table(i,1)  event_table(i+1,1) ];
            
        end
    end
    
    
    %== Merging event== %
    
    if( length( unique( raw_events_name ) ) < number_of_events ) 
       % Two or more events have the same name;
       events_index=containers.Map('KeyType','char','ValueType','any');
       for i = 1:(n-2)
            key=int2str(bi2de( event_table(i,2:end),'right-msb' ));
            key_name=events_names( key);
            next_key=events_names( int2str(bi2de( event_table(i+1,2:end),'right-msb' )));

            %if two consecutive events have the same name, we can merge them
            if(  strcmp(key_name, next_key) )
                event_table(i+1,1)=event_table(i,1);
            else  
                if( isKey (events_index,key) ) 
                    events_index(key)=[  events_index(key) ;  event_table(i,1) event_table(i+1,1)  ];
                else
                    events_index(key)=[ event_table(i,1)  event_table(i+1,1) ];
                end
            end
        end
       i=n-1;
       % adding the last event
       if( isKey (events_index,key) ) 
           events_index(key)=[  events_index(key) ;  event_table(i,1) event_table(i+1,1)  ];
       else
           events_names(key)=key;
           events_index(key)=[ event_table(i,1)  event_table(i+1,1) ];
       end
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
        duration_info = sprintf(', avg duration=%1.3f sec.', mean((newEvents(ievt).samples(2, :) - newEvents(ievt).samples(1, :)))/12.5);
        message{end+1} = sprintf(' - %s : %d trials. 1st trial at %1.3f sec.%s', ...
                                 newEvents(ievt).label, length(newEvents(ievt).epochs), ...
                                 newEvents(ievt).samples(1,1)/12.5, duration_info);
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
        newEvents(iNew).samples = round(newEvents(iNew).samples);
        newEvents(iNew).times   = newEvents(iNew).samples ./ sFile.prop.sfreq;
    else
        newEvents(iNew).samples = round(newEvents(iNew).times .* sFile.prop.sfreq);
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
        sFile.events(iEvt).samples    = [sFile.events(iEvt).samples, newEvents(iNew).samples];
        sFile.events(iEvt).epochs     = [sFile.events(iEvt).epochs, newEvents(iNew).epochs];
        sFile.events(iEvt).reactTimes = [sFile.events(iEvt).reactTimes, newEvents(iNew).reactTimes];
        % Sort by sample indices
        if (size(sFile.events(iEvt).samples, 2) > 1)
            [tmp__, iSort] = unique(sFile.events(iEvt).samples(1,:));
            sFile.events(iEvt).samples = sFile.events(iEvt).samples(:,iSort);
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
