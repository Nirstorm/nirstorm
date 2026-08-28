function varargout = process_nst_bootstrap_MEM( varargin )
% process_nst_bootstrap_MNE:  Compute source-localization using Bootstrap
% and MEM

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
%
% Copyright (c)2000-2018 University of Southern California & McGill University
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
% Authors: Edouard Delaire, 2022

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    % Description the process
    sProcess.Comment     = 'Boostrap MEM';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', '3D reconstruction'};
    sProcess.Index       = 1506;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'ressults'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;


    sProcess.options.label0.Comment = '<b>Bootstraps options:</b>';
    sProcess.options.label0.Type = 'label';

    sProcess.options.TimeSegmentNoise.Comment = 'Baseline Time window:';
    sProcess.options.TimeSegmentNoise.Type    = 'baseline';
    sProcess.options.TimeSegmentNoise.Value   = [];

    sProcess.options.timewindow_snr.Comment = 'Time window for SNR:';
    sProcess.options.timewindow_snr.Type    = 'poststim';
    sProcess.options.timewindow_snr.Value   =  [];

    sProcess.options.replacement.Comment = 'Do bootstrap with replacement ?';
    sProcess.options.replacement.Type    = 'checkbox';
    sProcess.options.replacement.Value   = 1;

    sProcess.options.combination.Comment = 'Combination value #:';
    sProcess.options.combination.Type    = 'value';
    sProcess.options.combination.Value   = {16, '', 0};

    sProcess.options.n_average.Comment = 'Number of average:';
    sProcess.options.n_average.Type    = 'value';
    sProcess.options.n_average.Value   = {50, '', 0};

    sProcess.options.target_quantile.Comment = {'median', '70 percentile', '80 percientile', '90 percentile', 'targeted SNR quantile:'; '0.5', '0.7', '0.8', '0.9',  ''} ;
    sProcess.options.target_quantile.Type    = 'radio_linelabel';
    sProcess.options.target_quantile.Value   = '0.5';

    % Time window for baseline
    sProcess.options.combination.Comment = 'Combination value #:';
    sProcess.options.combination.Type    = 'value';
    sProcess.options.combination.Value   = {16, '', 0};
    
    % Definition of the options for source reconstruction
    sProcess.options.label1.Comment = '<b>MEM options:</b>';
    sProcess.options.label1.Type = 'label';

    sProcess.options.thresh_dis2cortex.Comment = 'Reconstruction Field of view (distance to montage border)';
    sProcess.options.thresh_dis2cortex.Type    = 'value';
    sProcess.options.thresh_dis2cortex.Value   = {5, 'cm',2};
    sProcess.options.thresh_dis2cortex.Group   = 'MNE';


    sProcess.options.mem.Comment = {'panel_brainentropy', 'Source estimation options: '};
    sProcess.options.mem.Type    = 'editpref';
    sProcess.options.mem.Value   = be_main();

    sProcess.options.output_all.Comment = 'Save all bootstrap samples (warning: requires a lot of memory)';
    sProcess.options.output_all.Type    = 'checkbox';
    sProcess.options.output_all.Value   = 1;
    sProcess.options.output_all.Group   = 'output';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 

    OutputFiles = {};

    % Install/load brainentropy plugin
    [isInstalled, errMessage] = bst_plugin('Install', 'brainentropy', 1);
    if ~isInstalled
        bst_error(sprintf('Unable to install Brainentropy: %s', errMessage));
        return;
    end

    OutputFiles = process_nst_bootstrap_MNE('Run', sProcess, sInputs);
end
