function varargout = process_nst_import_evt_events( varargin )
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
    sProcess.Comment     = 'Import events from .evt table file';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS', 'Import'};
    sProcess.Index       = 1199;
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
    
    sProcess.options.label_cols.Comment = 'Events Names separated by commas';
    sProcess.options.label_cols.Type    = 'text';
    sProcess.options.label_cols.Value    = '';
    
    sProcess.options.last_event.Comment = 'Assume that last event has the same duration as the previous one';
    sProcess.options.last_event.Type = 'checkbox';
    sProcess.options.last_event.Value = 1; 
    
    sProcess.options.confirm_importation.Comment = 'Confirm importation';
    sProcess.options.confirm_importation.Type = 'checkbox';
    sProcess.options.confirm_importation.Value = 0; 
    sProcess.options.confirm_importation.Hidden  = 1;
    
    
    sProcess.options.preview.Comment =  {'process_nst_import_evt_events(''preview_importation'',iProcess,sfreq);' , '', 'Preview importation'} ;
    sProcess.options.preview.Type = 'button';
    sProcess.options.preview.Value   = [];
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    OutputFiles = {};
    
    isRaw = strcmpi(sInputs.FileType, 'raw');
    if ~isRaw
        bst_error('The process has to be run on raw file');
    end
    
    % Process inputs & do checks
    event_file  = sProcess.options.evtfile.Value{1};
    raw_events_name= nst_strsplit(sProcess.options.label_cols.Value, ','); 
    
    
    if ~exist(event_file, 'file')
        bst_error(['Event file not found: ' event_file]);
        return;
    end
    
    event_table = load(event_file);

    newEvents=compute(event_table,raw_events_name, sProcess.options.last_event.Value, sInputs(1));
    output_file = import_events(sProcess, sInputs(1), newEvents);

     if ~isempty(output_file)
         OutputFiles{end+1} = output_file;
     end
    
    
end

function newEvents =compute(event_table,raw_events_name,import_last, sInputs)

if nargin < 4
    sInputs = [];
end

    % we create a map for each events events_index(evt) will cointain every
    % time marker related to evt 
    
    events_names=containers.Map('KeyType','char','ValueType','char');
    events_index=containers.Map('KeyType','char','ValueType','any');
    
    n=length(event_table);
    number_of_events=length(raw_events_name);
    k=1;
    
    
    event_table= [ event_table(:,1) vect_b2d(event_table(:,2:end)) ];
    
    for i = 1:(n-1)

        key = int2str(event_table(i,2));
        
        if( isKey (events_index,key) ) 
            events_index(key)= [  events_index(key) ;  event_table(i,1)  event_table(i+1,1)  ];
        else
            if(~isempty(raw_events_name))
                if k > number_of_events 
                    bst_error('Not enough event name');
                    return;
                end
                events_names(key)=char(raw_events_name(k));
                
            else
                events_names(key)=num2str(k);
            end 
            
            k=k+1;
            events_index(key)= [ event_table(i,1)  event_table(i+1,1)];
            
        end
    end
    
    % importing the last event
    if import_last
        key = int2str(event_table(n,2));

        if( isKey(events_index,key) ) 
            evt=events_index(key);
            evt_duration=evt(end,2) - evt(end,1);
            events_index(key)= [  events_index(key) ;  event_table(end,1)  event_table(end,1)+evt_duration];
        else
            bst_report('Warning', 'process_nst_import_evt_events', sInputs, 'Not able to import last event');
        end
   end    
    
    
    %== Merging event== %
    
    if( length( unique( raw_events_name ) ) < number_of_events ) 
       % Two or more events have the same name;
       events_index=containers.Map('KeyType','char','ValueType','any');
       for i = 1:(n-2)
            key = int2str(event_table(i,2)); %(i,end:-1:2)
            key_name=events_names( key);
            next_key = events_names( int2str( event_table(i+1,2)));
            
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
        
       % adding the last event
       i=n-1;
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

time_offset=sFile.prop.times(1);

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
        newEvents(iNew).times   = time_offset+ newEvents(iNew).samples ./ sFile.prop.sfreq;
    else
        newEvents(iNew).samples = round(newEvents(iNew).times .* sFile.prop.sfreq);
        newEvents(iNew).times   = time_offset + newEvents(iNew).samples ./ sFile.prop.sfreq;
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


function y = b2d(x)
% Convert a binary array to a decimal number
% 
% Similar to bin2dec but works with arrays instead of strings and is found to be 
% rather faster
%
% From FileExchange: 
% https://www.mathworks.com/matlabcentral/fileexchange/26447-efficient-convertors-between-binary-and-decimal-numbers
%
% Copyright (c) 2010, Zacharias Voulgaris
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

z = 2.^(length(x)-1:-1:0);
y = sum(x.*z);

end

function y=vect_b2d(x)
    y=zeros(size(x,1),1);
    
    for i=1:size(x,1)
        y(i)=b2d(x(i,:));
    end
end

function preview_importation(iProcess,sfreq)
        % Get current process options
    global GlobalData;
    sProcess = GlobalData.Processes.Current(iProcess);
    
     file=sProcess.options.evtfile.Value;
     if ~isempty(file{1}) && exist(  file{1} ) 
        event_table = load(file{1});
        newEvents=compute(event_table,[],1);
  
        % Get existing specification figure
        hFig = findobj(0, 'Type', 'Figure', 'Tag', 'EventImport');
        % If the figure doesn't exist yet: create it
        if isempty(hFig)
        	hFig = figure(...
                'MenuBar',     'none', ...
                 'Toolbar',     'none', ...
                'NumberTitle', 'off', ...
                'Name',        sprintf('Preview of events importation'), ...
                'Tag',         'EventImport', ...
                'Units',       'Pixels');
        
        else
            clf(hFig);
            figure(hFig);
        end
        [path name ext]=fileparts(file{1});    
        message = ['<html> ' num2str( numel(newEvents))  ' events have been detected in  : ' name  ext  '<br />' ];
        for ievt=length(newEvents):-1:1
            tmp= sprintf(' - %s : %d trials. 1st trial at %1.3f sec , avg duration=%1.3f sec. <br />', ...
                                 newEvents(ievt).label, length(newEvents(ievt).epochs), ...
                                 newEvents(ievt).samples(1,1)/sfreq,  mean((newEvents(ievt).samples(2, :) - newEvents(ievt).samples(1, :)))/sfreq);
            message=[message tmp];                 
        end
        message=[message '</html>'];
        
        
        jLabel = javaObjectEDT('javax.swing.JLabel',message);
        [hcomponent,hcontainer] = javacomponent(jLabel,[0,0,600,80],hFig);




    end
        
end  
