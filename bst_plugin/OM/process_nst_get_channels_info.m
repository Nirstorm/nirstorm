function varargout = process_nst_get_channels_info( varargin )

% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2016 University of Southern California & McGill University
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
% Authors: Edouard Delaire, 2020; Thomas Vincent, 2015-2019


eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    sProcess.Comment     = 'Get channels info';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Optimal Montage'};
    sProcess.Index       = 1103;
    sProcess.Description = '';
    sProcess.isSeparator = 0; 
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    % Definition of the options
    % Description of the process
    sProcess.options.user_input.Comment = 'Output folder';
    sProcess.options.user_input.Type    = 'text';
    sProcess.options.user_input.Value    = '';

end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
    OutputFiles     = {sInputs.FileName};
    output_folder   = sProcess.options.user_input.Value;
    
    if isempty(output_folder)
        bst_error('Please enter a valid output folder. ');
        return
    end
    
    if length(unique({sInputs.SubjectName})) > 1
        bst_error('Cannot process multiple subjects in process')
        return
    end
    
    unique_input_conditions = unique({sInputs.Condition});
    
    for iCondition =  1:length(unique_input_conditions)
        iFile = find(strcmp({sInputs.Condition}, unique_input_conditions{iCondition}), 1);
        sChannel = in_bst_channel(sInputs(iFile).ChannelFile);
        
        all_S = [];
        all_D = [];
        
        SourcesLoc = NaN(500, 3);
        DetectorsLoc = NaN(500, 3);
        PairExists = false(500, 500);
        
        for iChannel = 1:length(sChannel.Channel)
            chName = sChannel.Channel(iChannel).Name;

            tokens = regexp(chName, 'S(\d+)D(\d+)', 'tokens', 'once');
            
            if ~isempty(tokens)
                s_idx = str2double(tokens{1});
                d_idx = str2double(tokens{2});
                
                all_S = [all_S, s_idx];
                all_D = [all_D, d_idx];
                
                if size(sChannel.Channel(iChannel).Loc, 2) >= 2
                    SourcesLoc(s_idx, :) = sChannel.Channel(iChannel).Loc(:,1)';
                    DetectorsLoc(d_idx, :) = sChannel.Channel(iChannel).Loc(:,2)';
                end
                
                PairExists(s_idx, d_idx) = true;
            end
        end
        
        if isempty(all_S)
            disp(['No channels with format SxDyWz found for condition : ', unique_input_conditions{iCondition}]);
            continue;
        end
        
        max_s = max(all_S);
        max_d = max(all_D);
        
        grid_text = cell(max_s + 1, max_d + 1);
        grid_class = cell(max_s + 1, max_d + 1);
        
        grid_text(:) = {'X'};
        grid_class(:) = {''};
        grid_text{1,1} = '';
        
        for d = 1:max_d
            grid_text{1, d+1} = sprintf('D%d', d);
        end
        for s = 1:max_s
            grid_text{s+1, 1} = sprintf('S%d', s);
        end
        
        for s = 1:max_s
            for d = 1:max_d
                if ~isnan(SourcesLoc(s,1)) && ~isnan(DetectorsLoc(d,1))
                    dist = euc_dist(SourcesLoc(s,:)', DetectorsLoc(d,:)');
                    dist_mm = dist * 1000; 
                    
                    grid_text{s+1, d+1} = sprintf('%.1f mm', dist_mm);
                    
                    if PairExists(s, d)
                        grid_class{s+1, d+1} = 'active';
                    end
                end
            end
        end
        
        safe_condition = regexprep(unique_input_conditions{iCondition}, '[^a-zA-Z0-9_]', '_');
        html_filename = fullfile(output_folder, sprintf('SD_matrix_%s.html', safe_condition));
        
        fid = fopen(html_filename, 'w');
        report_string = '';
        
        if fid == -1
            bst_error(['Impossible to create the file: ', html_filename]);
            return;
        end
        
        fprintf(fid, '<!DOCTYPE html>\n<html>\n<head>\n<style>\n');
        fprintf(fid, 'body { font-family: Arial, sans-serif; }\n');
        fprintf(fid, 'table { border-collapse: collapse; margin-top: 20px; }\n');
        fprintf(fid, 'th, td { border: 1px solid black; padding: 8px; text-align: center; }\n');
        fprintf(fid, 'th { background-color: #f2f2f2; font-weight: bold; }\n');
        fprintf(fid, '.active { background-color: lightgreen; font-weight: bold; }\n');
        fprintf(fid, '</style>\n</head>\n<body>\n');
        fprintf(fid, '<h3>Source-detector matrix (Get channels info process NIRSTORM) : </h3></h4>\n%s</h4>\n', safe_condition);
        
        fprintf(fid, '<table>\n');
        
        for row = 1:size(grid_text, 1)
            fprintf(fid, '  <tr>\n');
            for col = 1:size(grid_text, 2)
                val = grid_text{row, col};
                cls = grid_class{row, col};
                
                if row == 1 || col == 1
                    fprintf(fid, '    <th>%s</th>\n', val);
                else
                    if strcmp(cls, 'active')
                        fprintf(fid, '    <td class="active">%s</td>\n', val);
                    else
                        fprintf(fid, '    <td>%s</td>\n', val);
                    end
                end
            end
            fprintf(fid, '  </tr>\n');
        end
        
        fprintf(fid, '</table>\n</body>\n</html>');
        fclose(fid);

        
    end
end

function d = euc_dist(p1, p2)
    d = sqrt(sum((p1 - p2).^2));
end