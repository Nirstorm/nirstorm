function varargout = process_nst_deglitch( varargin )

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
% Authors: Thomas Vincent (2019)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>

    % Description the process
    sProcess.Comment     = 'Remove glitches';
    sProcess.Category    = 'Filter';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1403;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Remove-glitches';
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    
    sProcess.options.factor_std_grad.Comment = 'Variation threshold: ';
    sProcess.options.factor_std_grad.Type    = 'value';
    sProcess.options.factor_std_grad.Value   = {2.5, '* std(gradient)', 2}; % last if nb of subunit digits, use 0 for integer
    
%     sProcess.options.tag_only.Comment = 'Only mark glitches as events';
%     sProcess.options.tag_only.Type    = 'checkbox';
%     sProcess.options.tag_only.Value   =  0;

%TODO: apply only to NIRS channels
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function sInputs = Run(sProcess, sInputs) %#ok<DEFNU>
sInputs.A = Compute(sInputs.A, sProcess.options.factor_std_grad.Value{1});
end

%% ===== Compute =====
function [nirs_deglitched, glitch_flags] = Compute(nirs_sig, std_factor_thresh)

if nargin < 2
    std_factor_thresh = 2.5;
end

[nb_positions, nb_samples] = size(nirs_sig);
glitch_flags = detect_glitches(nirs_sig, std_factor_thresh);
nirs_deglitched = nirs_sig;
for ipos=1:nb_positions
    nirs_deglitched(ipos, glitch_flags(ipos, :)) = repmat(mean(nirs_sig(ipos, :)), ...
                                                          1, sum(glitch_flags(ipos, :)));
    if glitch_flags(ipos, end)
       nirs_deglitched(ipos, end) = mean(nirs_sig(ipos, nb_samples-2:nb_samples-1));
       glitch_flags(ipos, end) = 0;
    end
    if glitch_flags(ipos, 1)
       nirs_deglitched(ipos, 1) = mean(nirs_sig(ipos, 2:3));
       glitch_flags(ipos, 1) = 0;
    end
    % Maybe the following can be vectorized but there was something strange with
    % find function on 2D array
    glitch_i = find(glitch_flags(ipos, :));
    nirs_deglitched(ipos, glitch_i) = (nirs_sig(ipos, glitch_i-1) + nirs_sig(ipos, glitch_i+1)) ./ 2;
end
end

function glitch_flags = detect_glitches(nirs_sig, std_factor_thresh)

% nirs_sig: (nb_positions, nb_samples)

[nb_pos, nb_samples] = size(nirs_sig);

if nargin < 2
    std_factor_thresh = 2.5;
end

signal_tmp = [nirs_sig(:, 2)  nirs_sig  nirs_sig(:, end-1)]; % mirror edges
grad = diff(signal_tmp, 1, 2);
abs_grad = abs(grad);

% Robust standard deviation of absolute gradient
% -> clip values within 1% to 99% of total range, to remove extreme values
% then compute reference std of absolute gradient over these clipped values

sorted_signal = sort(signal_tmp, 2);
thresh_low = sorted_signal(:, round(nb_samples*0.01));
thresh_high = sorted_signal(:, round(nb_samples*0.99));
for ipos=1:size(signal_tmp, 1)
    signal_tmp(ipos, signal_tmp(ipos, :) < thresh_low(ipos)) = thresh_low(ipos);
    signal_tmp(ipos, signal_tmp(ipos, :) > thresh_high(ipos)) = thresh_high(ipos);
end
std_agrad = std(abs(diff(signal_tmp, 1, 2)), 0, 2);

glitch_canditates = [zeros(nb_pos, 1) (abs_grad > std_factor_thresh * std_agrad) .* sign(grad)];
glitch_flags = [abs(diff(glitch_canditates, 1, 2))==2 false(nb_pos, 1)] & glitch_canditates;
glitch_flags = glitch_flags(:, 2:end-1); % remove mirrored edges
end
