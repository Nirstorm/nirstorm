function varargout = process_nst_get_data_perform_2018( varargin )

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
% Authors: Thomas Vincent (2018)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Setup workshop PERFORM 2018';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS';
sProcess.Index       = 1902;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Workshop-PERFORM-Week-2018#download-data-sets-and-check-installation';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;

msg = 'Path to the downloaded and uncrompressed archive:';
sProcess.options.hint1.Comment = msg;
sProcess.options.hint1.Type = 'label';

% msg = '<I>in the given directory</I>';
% sProcess.options.hint2.Comment = msg;
% sProcess.options.hint2.Type = 'label';

SelectOptions = {...
    '', ...                            % Filename
    '', ...                            % FileFormat
    'save', ...                        % Dialog type: {open,save}
    'Select output folder...', ...     % Window title
    'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
    'single', ...                      % Selection mode: {single,multiple}
    'dirs', ...                        % Selection mode: {files,dirs,files_and_dirs}
    {{'.folder'}, '*.*'}, ... % Available file formats
    'MriOut'};                         % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,AnatIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}

% options.outputdir.Hidden = 1;
sProcess.options.inputdir.Comment = '';
sProcess.options.inputdir.Type    = 'filename';
sProcess.options.inputdir.Value   = SelectOptions;

sProcess.options.bst_dir.Hidden = 1;
sProcess.options.bst_dir.Comment = 'Bst home directory';
sProcess.options.bst_dir.Type    = 'text';
sProcess.options.bst_dir.Value   = bst_get('BrainstormUserDir');

sProcess.options.confirm_download.Hidden = 1;
sProcess.options.confirm_download.Comment = 'Confirm download';
sProcess.options.confirm_download.Type    = 'checkbox';
sProcess.options.confirm_download.Value   = 0;

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
global GlobalData;

OutputFiles = {};

data_dir = sProcess.options.inputdir.Value{1};
if ~exist(data_dir, 'dir')
    bst_error(sprintf('Input data folder not found: %s', data_dir));
end
bst_dir = sProcess.options.bst_dir.Value;

%% Create protocol if not already there
protocol_name = 'nirstorm_perform_2018';
% i_protocol = ;
if isempty(bst_get('Protocol', protocol_name))
    % Create new protocol
    gui_brainstorm('CreateProtocol', protocol_name, 0, 0);
end
gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', protocol_name));

bst_report('Start');

%% Check sample data sets
bst_progress('start', 'Data check','Checking package data...');
data_fns = get_data_file_names(data_dir);

if files_exist(data_fns)
    bst_report('Info', sProcess, sInputs, 'Input sample data ok');
else
    % TODO: show download URL
    bst_error(sprintf('Sample data not found in %s', data_dir));
    return;
end
bst_progress('stop');

%% Install template
bst_progress('start', 'Install template', 'Installing template Colin27 for NIRS...');
src_fn = fullfile(data_dir, 'Colin27_4NIRS.zip');
dest_fn = get_colin27_installed_fn(bst_dir);
copyfile(src_fn, dest_fn);
bst_progress('stop');

if exist(dest_fn, 'file')
    bst_report('Info', sProcess, sInputs, 'Template data installation ok');
else
    bst_error('Template data installation failed');
end

%% Install fluences
% Uncompress archive
bst_progress('start', 'Unzip fluences', 'Unzipping fluence files...');
fluence_zip_fn = data_fns{end};
[root, fluence_base_folder, ext] = fileparts(fluence_zip_fn);
fluence_folder = fullfile(data_dir, fluence_base_folder);
unzip(fluence_zip_fn, fluence_folder);
bst_progress('stop');

% Import fluences from uncompressed archive
fluence_fns = get_OM_fluence_fns(bst_dir);
fluence_install_dir = fileparts(fluence_fns{1});
if ~exist(fluence_install_dir, 'dir')
    mkdir(fluence_install_dir);
end
bst_progress('start', 'Install fluences', 'Installing fluence files...', 1, length(fluence_fns));
for ifluence=1:length(fluence_fns)
    dest_fn = fluence_fns{ifluence};
    [root, bfn, ext] = fileparts(dest_fn);
    src_fn = fullfile(fluence_folder, [bfn ext]);
    copyfile(src_fn, dest_fn);
    bst_progress('inc', 1);
end
rmdir(fluence_folder, 's');
bst_progress('stop');

if files_exist(fluence_fns)
    bst_report('Info', sProcess, sInputs, 'Fluence data installation ok');
else
    bst_error('Fluence data installation failed');
end

bst_progress('start', 'Check dependencies', 'Checking dependencies...', 1, 3);
%% Check CPLEX
try
    cplx = Cplex();
    cplex_version = strsplit(cplx.getVersion(), '.');
    if str2double(cplex_version(1)) < 12 || ...
            (length(cplex_version) > 1 && str2double(cplex_version(1)) < 3)
        bst_error(['CPLEX >12.3 required. See ' cplex_url]);
        return
    end
    bst_report('Info', sProcess, sInputs, 'Cplex >= v12.3 found');
catch
    bst_report('Error', sProcess, sInputs, 'Cplex >= v12.3 not found.<BR>Required for optimal montage computation.');
end
bst_progress('inc', 1);

if ~license('test', 'Curve_Fitting_Toolbox')
    bst_report('Warning', sProcess, sInputs, 'Curve Fitting Toolbox not available. Motion correction will not work.');
else
    if isempty(which('csaps'))
        bst_report('Warning', sProcess, sInputs, ...
                   ['Curve Fitting Toolbox OK but function csaps not found.<BR>' ...
                    'Try refreshing matlab cache using command: rehash toolboxcache']);
    else
        bst_report('Info', sProcess, sInputs, 'Curve Fitting Toolbox found');
    end
end
bst_progress('inc', 1);

bst_progress('stop');

bst_interactive = ~isempty(GlobalData) && isfield(GlobalData, 'Program') && ...
                  ~isempty(GlobalData.Program) && ...
                  (~isfield(GlobalData.Program, 'isServer') || ...
                   ~GlobalData.Program.isServer);
if bst_interactive
    bst_report('Open', 'current');
end

end

function data_fns = get_data_file_names(data_dir)
% workshop_dir = 'Nirstorm_workshop_PERFORM_2018';
if nargin < 1
    data_dir = '';
end
data_fn_parts = {{data_dir, 'channel_space_analysis','sample_nirs', 'data', 'fiducials.txt'}, ...
             {data_dir, 'channel_space_analysis','sample_nirs', 'data', 'optodes.txt'}, ...
             {data_dir, 'channel_space_analysis','sample_nirs', 'data', 'S01_Block_FO_LH_Run01.nirs'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'nobias_JegAud.nii'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'JegAud_head.gii'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'JegAud_Lhemi.gii'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'JegAud_Rhemi.gii'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'JegAud_Lwhite.gii'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'JegAud_Rwhite.gii'},...
             {data_dir, 'channel_space_analysis','sample_nirs', 'anatomy', 'JegAud.APC'},...
             {data_dir, 'cortical_reconstruction','sample_NIRS_Reconstruction.zip'},...
             {data_dir, 'Colin27_4NIRS.zip'}, {data_dir, 'fluences_for_OM.zip'}};
data_fns = cellfun(@(pt) fullfile(pt{:}), data_fn_parts, 'UniformOutput', false);
end

function vertex_ids = get_OM_scout_vertex_ids()
vertex_ids = [41 147 150 156 157 160 161 568 572 576 592 594 598 606 608 609 610 611 612 616 618 619 621 624 625 626 627 632 633 636 637 638 2233 2239 2257 2281 2287 2294 2298 2302 2317 2321 2325 2329 2333 2336 2339 2345 2347 2367 2371 2377 2381 2385 2391 2393 2395 2407 2410 2416 2417 2418 2419 2420 2421 2430 2432 2434 2438 2439 2440 2441 2447 2449 2450 2451 2452 2453 2454 2455 2456 2465 2469 2471 2472 2473 2474 2479 2481 2482 2483 2484 2485 2486 2487 2488 2495 2499 2500 2501 2502 2504 2506 2507 2508 2509 2510 2511 2512 2513 2518 2522 2523 2524 2526 2528 2529 2530 2531 2532 2533 2538 2539 2544 2545 2546 2547 2552 8939 8979 8991 9009 9021 9027 9061 9064 9069 9077 9083 9100 9109 9117 9129 9161 9165 9167 9179 9183 9189 9197 9201 9213 9215 9217 9225 9243 9247 9252 9259 9265 9267 9271 9273 9282 9287 9294 9300 9304 9306 9314 9348 9352 9358 9364 9366 9371 9378 9380 9382 9386 9390 9394 9400 9404 9410 9414 9433 9437 9439 9445 9447 9449 9453 9459 9463 9465 9467 9477 9481 9485 9489 9491 9496 9532 9534 9536 9542 9544 9546 9550 9552 9558 9560 9564 9568 9576 9580 9582 9584 9586 9592 9594 9596 9616 9620 9622 9626 9660 9662 9664 9670 9674 9677 9678 9679 9680 9681 9682 9693 9697 9699 9703 9705 9706 9707 9708 9709 9710 9711 9712 9713 9714 9715 9716 9717 9718 9719   9720   9735   9739   9745   9747   9749   9750   9751   9752   9753   9754   9755   9756   9763   9765   9767   9771   9773   9774   9775   9776   9777   9778   9779   9780   9781   9782   9783   9784   9785   9786   9787   9788   9789   9803   9807   9813   9816   9818   9819   9820   9821   9822   9823   9824   9825   9831   9833   9835   9841   9842   9843   9844   9845   9846   9847   9848   9849   9850   9851   9852   9853   9854   9855   9856   9857   9867   9871   9875   9879   9880   9881   9882   9883   9884   9885   9886   9891   9893   9895   9901   9902   9903   9904   9905   9906   9907   9908   9909   9910   9911   9912   9913   9914   9915   9916   9927   9931   9933   9937   9938   9939   9940   9941   9942   9943   9944   9948   9950   9958   9959   9960   9961   9962   9963   9964   9965   9966   9967   9968   9969   9970   9971   9972   9973   9980   9982   9984   9988   9989   9990   9991   9992   9993   9994   9995   9998  10002  10010  10011  10012  10013  10014  10015  10016  10017  10018  10019  10020  10021  10022  10023  10024  10030  10032  10034  10038  10039  10040  10041  10042  10043  10044  10049  10053  10054  10055  10056  10057  10058  10059  10060  10061  10062  10063  10064  10065  10066  10073  10075  10079  10080  10081  10082  10083  10084  10085  10090  10094  10095  10096  10097  10098  10099  10100  10101  10102  10103  10104  10105  10106  10112  10114  10118  10119  10120  10121  10122  10123  10127  10129  10130  10131  10132  10133  10134  10135  10136 10137  10138  10145  10151  10152 10153 10154 10159 10161 10162 10163 10164 10165 10166 10167 10168 10175 10179 10180 10181 10186 10187 10188 10189 10190 10191 10202 10203 10208 10209];
end

function fluence_fns = get_OM_fluence_fns(bst_dir)
vertex_ids = get_OM_scout_vertex_ids();
local_repository = fullfile(bst_dir, 'defaults', 'nirstorm', 'fluence', 'MRI__Colin27_4NIRS');
fluence_bfns = arrayfun(@(vid) process_nst_import_head_model('get_fluence_fn', vid, 685),...
                        vertex_ids, 'UniformOutput', false);
fluence_fns = cellfun(@(bfn) fullfile(local_repository, bfn), fluence_bfns, 'UniformOutput', false);
end

function template_fn = get_colin27_installed_fn(bst_dir)
template_fn = fullfile(bst_dir, 'defaults', 'anatomy', 'Colin27_4NIRS.zip');
end

function flag = files_exist(fns)
flag = all(cellfun(@(fn) exist(fn, 'file'), fns));
end
