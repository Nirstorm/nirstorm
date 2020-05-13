function varargout = process_nst_extract_sensitivity_from_head_model( varargin )

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
% Authors: Thomas Vincent, ZhengChen Cai (2017)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Extract sensitivity surfaces from head model';
    sProcess.FileTag     = '';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'NIRS';
    sProcess.Index       = 1405;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results', 'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

% Save the new head model
sStudy = bst_get('Study', sInputs.iStudy);

if isempty(sStudy.iHeadModel)
    bst_error('No head model found. Consider process "Compute head model from fluence"');
    return;
end

head_model = in_bst_headmodel(sStudy.HeadModel(sStudy.iHeadModel).FileName);
ChannelMat = in_bst_channel(sInputs(1).ChannelFile);

if ~strcmp(head_model.HeadModelType, 'surface')
    bst_error('Extraction only works for surface head model');
    return;
end

if ndims(head_model.Gain) ~= 3
    % TODO: better test shape consistency
   bst_error('Bad shape of gain matrix, must be nb_pairs x nb_wavelengths x nb_vertices');
   return;
end
 
montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel);

src_coords = montage_info.src_pos;
det_coords = montage_info.det_pos;

nb_sources = size(src_coords, 1);
nb_dets = size(det_coords, 1);
sensitivity_surf = head_model.Gain;
pair_names = head_model.pair_names;
nb_nodes = size(sensitivity_surf, 3);
% Save sensitivities
for iwl=1:size(sensitivity_surf, 2)
    sensitivity_surf_sum = sum(sensitivity_surf(:, iwl, :),  1) ;
    [sStudy, ResultFile] = add_surf_data(repmat(squeeze(sensitivity_surf_sum), [1,2]), [0 1], ...
                                         head_model, ['Summed sensitivities - WL' num2str(iwl)], ...
                                         sInputs.iStudy, sStudy,  ...
                                        'sensitivity imported from MCXlab');
    
    OutputFiles{end+1} = ResultFile;

end

if nb_dets < 100
    time = 1:(nb_sources*100 + nb_dets);
    for iwl=1:size(sensitivity_surf, 2)
        %sens_tmp = zeros(nb_nodes, length(time)) - 1;
        sens_tmp = zeros(nb_nodes, length(time));
        for ipair=1:size(sensitivity_surf, 1)
            [src_id, det_id] = nst_unformat_channel([pair_names{ipair} 'WL0']);
            sens_tmp(:, det_id + src_id*100) = squeeze(sensitivity_surf(ipair, iwl, :)); %source id will be minutes, det_it will be seconds
        end
        [sStudy, ResultFile] = add_surf_data(sens_tmp, time, ...
            head_model, ['Sensitivities - WL' num2str(iwl)], ...
            sInputs.iStudy, sStudy, 'Sensitivity import from template'); %TODO better denomitation
        OutputFiles{end+1} = ResultFile;
    end
end

% for iwl=1:size(sensitivity_surf, 2)
%     for ipair=1:size(sensitivity_surf, 1)
%         [sStudy, ResultFile] = add_surf_data(repmat(squeeze(sensitivity_surf(ipair, iwl, :)), [1,2]), [0 1], ...
%                                              newHeadModel, ['Sensitivity - WL' num2str(iwl) ' ' pair_names{ipair}], ...
%                                              iStudy, sStudy, 'pair sensitivity imported from MCXlab');
%         OutputFiles{end+1} = ResultFile;
%     end
% end

end


function [sStudy, ResultFile] = add_surf_data(data, time, head_model, name, ...
                                              iStudy, sStudy, history_comment)
                                          
%% Save a cortical map to brainstorm with given data

    ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                            ['results_' protect_fn_str(name)]);
                        
    % ===== CREATE FILE STRUCTURE =====
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    ResultsMat.ImageGridAmp = data;
    ResultsMat.Time          = time;
    ResultsMat.DataFile      = [];
    ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName; 
    ResultsMat.HeadModelType = head_model.HeadModelType;
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = file_short(head_model.SurfaceFile);
    ResultsMat.GridLoc    = [];
    ResultsMat.GridOrient = [];
    ResultsMat.nAvg      = 1;
    % History
    ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
    % Save new file structure
    bst_save(ResultFile, ResultsMat, 'v6');
    % ===== REGISTER NEW FILE =====
    % Create new results structure
    newResult = db_template('results');
    newResult.Comment       = name;
    newResult.FileName      = file_short(ResultFile);
    newResult.DataFile      = ''; %sInputs.FileName;
    newResult.isLink        = 0;
    newResult.HeadModelType = ResultsMat.HeadModelType;
    % Add new entry to the database
    iResult = length(sStudy.Result) + 1;
    sStudy.Result(iResult) = newResult;
    % Update Brainstorm database
    bst_set('Study', iStudy, sStudy);
end


function sfn = protect_fn_str(s)
sfn = strrep(s, ' ', '_');
end
