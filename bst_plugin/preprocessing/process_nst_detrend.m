function varargout = process_nst_detrend( varargin )

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
% Authors: Delaire Edouard (2022)
%
%
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Remove slow fluctuations';
sProcess.FileTag     = @GetFileTag; 
sProcess.Category    = 'Filter';
sProcess.SubGroup    = {'NIRS', 'Pre-process'};
sProcess.Index       = 1309;
sProcess.isSeparator = 0;
sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/Optode-separations';

% Definition of the input accepted by this process
sProcess.InputTypes  = {'data','raw'};
sProcess.OutputTypes = {'data','data'};


sProcess.options.filter_model.Comment = {'Polynome', 'DCT', 'Trend modelisation: '; ...
                                        'legendre', 'DCT',''};
sProcess.options.filter_model.Type    = 'radio_linelabel';
sProcess.options.filter_model.Value   = 'DCT';
sProcess.options.filter_model.Controller = struct('legendre','legendre','DCT','DCT' );

    
sProcess.options.option_period.Comment = 'Miminum Period (0= only linear detrend):';
sProcess.options.option_period.Type    = 'value';
sProcess.options.option_period.Value   = {200, 's', 0};
sProcess.options.option_period.Class = 'DCT';
    

sProcess.options.poly_order.Comment = 'Polynome order:';
sProcess.options.poly_order.Type    = 'value';
sProcess.options.poly_order.Value   = {3, '', 0};
sProcess.options.poly_order.Class = 'legendre';


sProcess.options.option_keep_mean.Comment = 'Keep the mean';
sProcess.options.option_keep_mean.Type    = 'checkbox';
sProcess.options.option_keep_mean.Value   = 0;
    
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 1;
sProcess.nOutputs    = 1;
end

%% ===== FORMAT COMMENT =====
function [Comment, fileTag] = FormatComment(sProcess)
    if strcmp(sProcess.options.filter_model.Type, 'legendre')
        Comment = 'Detend (polynoial)';
    elseif strcmp(sProcess.options.filter_model.Type, 'DCT')
        Comment = 'Detend (DCT)';
    else
        Comment = 'Detend';
    end

    fileTag = 'detrend';
end
%% ===== GET FILE TAG =====
function fileTag = GetFileTag(sProcess)
    [Comment, fileTag] = FormatComment(sProcess);
end

function sInput = Run(sProcess, sInput)

    keep_mean = sProcess.options.option_keep_mean.Value;
    ChannelMat = in_bst_channel(sInput.ChannelFile);
    nirs_ichans = sInput.ChannelFlag ~= -1 & strcmpi({ChannelMat.Channel.Type}, 'NIRS')';
    
    
    Y =  sInput.A(nirs_ichans,:)';
    
    model = nst_glm_initialize_model(sInput.TimeVector);
    model = nst_glm_add_regressors(model,'constant');
    
    switch sProcess.options.filter_model.Value
        case 'DCT'
            model = nst_glm_add_regressors(model,'linear');  
            if sProcess.options.option_period.Value{1} > 0
               period   = sProcess.options.option_period.Value{1};
               model    = nst_glm_add_regressors(model,'DCT',[1/sInput.TimeVector(end) 1/period],{'LFO'});  
            end
        case 'legendre'
            model = nst_glm_add_regressors(model,'legendre', sProcess.options.poly_order.Value{1} );  
    end
    
    
    [B,proj_X] = nst_glm_fit_B(model,Y, 'SVD');
    
    if keep_mean
        Y = Y - model.X(:,2:end)*B(2:end,:);
    else    
        Y = Y - model.X*B;
    end
    
    
    sInput.A(nirs_ichans,:)     = Y';
    sInput.CommentTag           = FormatComment(sProcess);

    History                     = 'Remove Linear trend';
    if strcmp(sProcess.options.filter_model.Value,'DCT') && sProcess.options.option_period.Value{1}
        History                 = [History sprintf('(using DCT of period > %ds)',period)];
      
    elseif strcmp(sProcess.options.filter_model.Value,'legendre') && sProcess.options.option_period.Value{1}
        History                 = [History sprintf('(using %d order polynomial)',sProcess.options.poly_order.Value{1})];
    end    
    sInput.HistoryComment = History;

end    