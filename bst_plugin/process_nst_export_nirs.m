function varargout = process_nst_export_nirs( varargin )

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
% Authors: Thomas Vincent (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Export to .nirs (HOMer format)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1007;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % Definition of the options    

    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'save', ...                        % Dialog type: {open,save}
        'Select output folder...', ...     % Window title
        'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'dirs', ...                        % Selection mode: {files,dirs,files_and_dirs}
        {{'.folder'}, '*.*'}, ...          % Available file formats
        'ImageOut'};                         % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    % Option definition
    sProcess.options.outputdir.Comment = 'Output folder for nirs file';
    sProcess.options.outputdir.Type    = 'filename';
    sProcess.options.outputdir.Value   = SelectOptions;
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

for iInput=1:length(sInputs)
    % Load channel file
    ChanneMat = in_bst_channel(sInputs(iInput).ChannelFile);
    
    % Load recordings
    if strcmp(sInputs(iInput).FileType, 'data')     % Imported data structure
        sDataIn = in_bst_data(sInputs(iInput).FileName);
        events = sDataIn.Events;
        time = sDataIn.Time';
    elseif strcmp(sInputs(iInput).FileType, 'raw')  % Continuous data file
        sDataIn = in_bst(sInputs(iInput).FileName, [], 1, 1, 'no');
        sDataRaw = in_bst_data(sInputs(iInput).FileName, 'F');
        events = sDataRaw.F.events;
        time = sDataIn.Time';
    end
    
    nirs = Compute(ChanneMat, sDataIn.F, time, events);
    
    cond_name = sInputs(iInput).Condition;
    if strcmp(cond_name(1:4), '@raw')
        cond_name = cond_name(5:end);
    end
    
    if strfind(sInputs(iInput).Comment, 'Link to raw file')
        suffix = '_raw';
    else
        suffix = ['_' sInputs(iInput).Comment];
    end
    if ~exist(sProcess.options.outputdir.Value{1}, 'dir')
        bst_error(['Output folder "' sProcess.options.outputdir.Value{1} '" does not exist']);
    end
    nirs_fn = fullfile(sProcess.options.outputdir.Value{1}, ...
                       protect_fn_str([cond_name suffix '.nirs']) );
    save(nirs_fn, '-struct', 'nirs');
end
end


function nirs = Compute(channel_def, data, time, events)
if size(time, 1) == 1;
    nirs.t = time';
else
    nirs.t = time;
end
if size(data, 2) ~= length(channel_def.Channel)
    data = data';
end

montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel);
pair_ichans = montage_info.pair_ichans;
src_coords = montage_info.src_pos;
det_coords = montage_info.det_pos;
pair_sd_indexes =  montage_info.pair_sd_indexes;

if isfield(channel_def.Nirs, 'Wavelengths')
    nirs.SD.Lambda = channel_def.Nirs.Wavelengths;
else
    nirs.SD.Lambda = channel_def.Nirs.Hb;
end
ichan = 1;
for imeasure=1:length(nirs.SD.Lambda)
    for ipair=1:size(pair_sd_indexes, 1)
        nirs.ml(ichan, 1) = pair_sd_indexes(ipair, 1);
        nirs.ml(ichan, 2) = pair_sd_indexes(ipair, 2);
        nirs.ml(ichan, 3) = 1;
        nirs.ml(ichan, 4) = imeasure;
        nirs.d(:, ichan) = data(:, pair_ichans(ipair, imeasure));
        ichan = ichan + 1;
    end
end

nirs.SD.MeasList = nirs.ml;
nirs.SD.SrcPos = src_coords;
nirs.SD.nSrcs = size(src_coords, 1);
nirs.SD.DetPos = det_coords;
nirs.SD.nDets = size(det_coords, 1);
nirs.aux = data(:, strcmp({channel_def.Channel.Type}, 'NIRS_AUX'));
nirs.events = events;
end



function fn = protect_fn_str(sfn)
fn = strrep(sfn, ' | ', '_');
fn = strrep(fn, ' ', '_');
fn = strrep(fn, '"', '');
fn = strrep(fn, ':', '_');
fn = strrep(fn, '(', '_');
fn = strrep(fn, ')', '_');
fn = strrep(fn, '[', '_');
fn = strrep(fn, ']', '_');
fn = strrep(fn, '!', '');
fn = strrep(fn, '__', '_');
end
