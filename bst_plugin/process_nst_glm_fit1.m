function varargout = process_nst_glm_fit1( varargin )
% process_compute_glm: compute the glm : find B such as Y = XB +e with X
% 
% OlS_fit use an ordinary least square algorithm to find B : B= ( X^{T}X)^{-1} X^{T} Y 
% AR-IRLS : Details about AR-IRLS algorithm can be found here :
%   http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3756568/ 
% 
% 
% Further update : use more sophisticated method to fit the B
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
% Authors: Edouard Delaire, Thomas Vincent 2018-2019
%
% TODO: use different output types for channel-space and surface analyses
%       -> more convenient for contrast computation afterwards, to be able to map
%       input to output types
%
eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'GLM - 1st level design and fit';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1601;
    sProcess.isSeparator = 0;

    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw', 'results'};
    sProcess.OutputTypes = {'data', 'data', 'results'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.label0.Comment = '<U><B>Signal Information</B></U>:';
    sProcess.options.label0.Type    = 'label';
    
    sProcess.options.hpf_low_cutoff.Comment = 'Applied High-pass filter ( 0 if no high-pass filter have been applied): ';
    sProcess.options.hpf_low_cutoff.Type    = 'value';
    sProcess.options.hpf_low_cutoff.Value   = {0.01, 'Hz', 2};
    
    sProcess.options.trim_start.Comment = 'Ignore starting signal: ';
    sProcess.options.trim_start.Type    = 'value';
    sProcess.options.trim_start.Value   = {0, 'sec', 2};

    sProcess.options.label3.Comment = '<U><B>Design Matrix</B></U>:';
    sProcess.options.label3.Type    = 'label';
    
    sProcess.options.stim_events.Comment = 'Stimulation events: ';
    sProcess.options.stim_events.Type    = 'text';
    sProcess.options.stim_events.Value   = '';
    
    sProcess.options.hrf_model.Comment = 'HRF model: ';
    sProcess.options.hrf_model.Type    = 'combobox';
    sProcess.options.hrf_model.Value   = {1, fieldnames(get_hrf_types())};
    
    sProcess.options.trend.Comment = 'Add constant regressor ';
    sProcess.options.trend.Type    = 'checkbox';
    sProcess.options.trend.Value   =  1;
    
    % Separator
    sProcess.options.separator00.Type = 'separator';
    sProcess.options.separator00.Comment = ' ';

    sProcess.options.label1.Comment = '<U><B>Optimization Method</B></U>:';
    sProcess.options.label1.Type    = 'label';
    sProcess.options.label1.Hidden   = 1;
    
    sProcess.options.fitting.Type    = 'radio_line';
    sProcess.options.fitting.Comment   = {'OLS', 'IRLS(not implemented)','' };
    sProcess.options.fitting.Value   = 1;
    sProcess.options.fitting.Hidden   = 1;
    
    sProcess.options.label2.Comment = '<U><B>Serial Correlation Preprocessing</B></U>:';
    sProcess.options.label2.Type    = 'label';
    
    sProcess.options.statistical_processing.Type    = 'radio_line';
    sProcess.options.statistical_processing.Comment   = {'Pre-coloring', 'Pre-whitenning','Method : '};
    sProcess.options.statistical_processing.Value   = 1;
    
    sProcess.options.output_cmt0.Comment = '<B>Pre-coloring Options</B>:';
    sProcess.options.output_cmt0.Type    = 'label';
    sProcess.options.output_cmt0.Hidden   = 1;
    
    sProcess.options.output_cmt1.Comment = '<B>Pre-whitenning Options</B>:';
    sProcess.options.output_cmt1.Type    = 'label';
    
    sProcess.options.noise_model.Type    = 'radio_line';
    sProcess.options.noise_model.Comment   = {'AR(1)', 'AR(p)','Model of the noise : '};
    sProcess.options.noise_model.Value   = 1;
    
    sProcess.options.separator1.Type = 'separator';
    sProcess.options.separator1.Comment = ' ';

    sProcess.options.output_cmt.Comment = '<U><B>Extra outputs</B></U>:';
    sProcess.options.output_cmt.Type    = 'label';
        
    sProcess.options.save_betas.Comment = 'Beta maps';
    sProcess.options.save_betas.Type    = 'checkbox';
    sProcess.options.save_betas.Value   =  0;

    sProcess.options.save_residuals.Comment = 'Residuals';
    sProcess.options.save_residuals.Type    = 'checkbox';
    sProcess.options.save_residuals.Value   =  0;
    
    sProcess.options.save_fit.Comment = 'Fit';
    sProcess.options.save_fit.Type    = 'checkbox';
    sProcess.options.save_fit.Value   =  0;
end



%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
%     fitting_choice=cell2mat(sProcess.options.fitting.Value(1));
%     if( fitting_choice == 1 )
%         Comment=[ Comment ' OLS fit'];
%     elseif( fitting_choice == 2)
%         Comment=[ Comment ' AR-IRLS fit'];
%     end    
end

function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    
    OutputFiles = bst_process('CallProcess', 'process_nst_glm_fit', sInput,   db_template('data'),   sProcess.options);

end
