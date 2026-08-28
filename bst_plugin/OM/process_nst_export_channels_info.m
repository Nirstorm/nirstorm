function varargout = process_nst_export_channels_info( varargin )
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
% Authors: Jean-Eudes Bornert, 2026 ; Edouard Delaire, 2026
eval(macro_method);
end

%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() 
    sProcess.Comment     = 'Export channels info';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Optimal Montage'};
    sProcess.Index       = 1110;
    sProcess.Description = '';
    sProcess.isSeparator = 0; 
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Definition of the options
    sProcess.options.include_images.Comment       = 'Include image(s) to the report';
    sProcess.options.include_images.Type          = 'checkbox';
    sProcess.options.include_images.Controller    = 'include_image';
    sProcess.options.include_images.Value         = 1;
    
    % JSON sidecar file
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'open', ...                           % Dialog type: {open,save}
        'Import image...', ...                % Window title
        'ExportImage', ...                    % LastUsedDir
        'single', ...                         % Selection mode
        'files', ...                          % Selection mode
        {{'.png'}, {'PNG document (*.png)'}, 'PNG'}, ...
        ''};
    
    sProcess.options.image1.Type    = 'filename';
    sProcess.options.image1.Comment = 'Image 1:';
    sProcess.options.image1.Value   = SelectOptions;
    sProcess.options.image1.Class   = 'include_image';
    
    sProcess.options.image2.Type    = 'filename';
    sProcess.options.image2.Comment = 'Image 2:';
    sProcess.options.image2.Value   = SelectOptions;
    sProcess.options.image2.Class   = 'include_image';
    
    sProcess.options.save_to_file.Comment       = 'Save report to file';
    sProcess.options.save_to_file.Type          = 'checkbox';
    sProcess.options.save_to_file.Controller    = 'is_save';
    sProcess.options.save_to_file.Value         = 1;
    sProcess.options.save_to_file.Group         = 'output';
    
    % JSON sidecar file
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'save', ...                           % Dialog type: {open,save}
        'Export PDF report...', ...           % Window title
        'ExportImage', ...                    % LastUsedDir
        'single', ...                         % Selection mode
        'files', ...                          % Selection mode
        {{'.pdf'}, {'PDF document (*.pdf)'}, 'PDF'}, ...
        ''};
    
    sProcess.options.user_input.Type    = 'filename';
    sProcess.options.user_input.Comment = 'Report file';
    sProcess.options.user_input.Value   = SelectOptions;
    sProcess.options.user_input.Group   = 'output';
    sProcess.options.user_input.Class   = 'is_save';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) 
    OutputFiles     = {sInputs.FileName};
    pdf_filename    = sProcess.options.user_input.Value{1};
    isSaveToFile    = sProcess.options.save_to_file.Value;
    
    if isSaveToFile && isempty(pdf_filename)
        bst_error('Please enter a valid output file.');
        return;
    end
    
    image_path_1 = sProcess.options.image1.Value{1}; 
    image_path_2 = sProcess.options.image2.Value{1}; 
    isIncludeImage = sProcess.options.include_images.Value;
    
    if isIncludeImage && ((~isempty(image_path_1) && ~file_exist(image_path_1)) || (~isempty(image_path_2) && ~file_exist(image_path_2)))
        bst_error('Please enter a valid image file.');
        return;
    end
    
    if length(unique({sInputs.SubjectName})) > 1
        bst_error('Cannot process multiple subjects in process');
        return;
    end
    
    unique_input_conditions = unique({sInputs.Condition});
    
    for iCondition = 1:length(unique_input_conditions)
        iFile = find(strcmp({sInputs.Condition}, unique_input_conditions{iCondition}), 1);
        sChannel = in_bst_channel(sInputs(iFile).ChannelFile);
        
        all_S = [];
        all_D = [];
        
        SourcesLoc = NaN(500, 3); SourceEEGName = cell(1, 500);
        DetectorsLoc = NaN(500, 3); DetectorEEGName = cell(1, 500);
        PairExists = false(500, 500);
        
        for iChannel = 1:length(sChannel.Channel)
            chName = sChannel.Channel(iChannel).Name;
            chComment = sChannel.Channel(iChannel).Comment;
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
                if ~isempty(chComment) && contains(chComment, '=>')
                    tmp = strsplit(chComment, '=>');
                    SourceEEGName{s_idx} = [' - ', strrep(tmp{1}, ' ', '')];
                    DetectorEEGName{d_idx} = [' - ',  strrep(tmp{2}, ' ', '')];
                end
            end
        end
        
        if isempty(all_S)
            disp(['No channels with format SxDyWz found for condition : ', unique_input_conditions{iCondition}]);
            continue;
        end
        
        title_montage = strrep(regexprep(unique_input_conditions{iCondition}, '[^a-zA-Z0-9_]', '_'), '_', ' ');    
        max_s = max(all_S);
        max_d = max(all_D);
        
        grid_text = cell(max_s + 1, max_d + 1);
        grid_class = cell(max_s + 1, max_d + 1);
        
        grid_text(:) = {'X'};
        grid_class(:) = {''};
        grid_text{1,1} = '';
        
        for d = 1:max_d
            grid_text{1, d+1} = sprintf('D%d%s', d, DetectorEEGName{d});
        end
        
        for s = 1:max_s
            grid_text{s+1, 1} = sprintf('S%d%s', s, SourceEEGName{s});
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
        
        report_string = '';
        report_string = [ report_string, sprintf('<!DOCTYPE html>\n<html>\n<head>\n<meta charset="utf-8">\n<style>\n')];
        report_string = [ report_string, sprintf('@page { size: A4 landscape; margin: 8mm; }\n')];
        report_string = [ report_string, sprintf('* { box-sizing: border-box; }\n')];
        report_string = [ report_string, sprintf('body { font-family: Arial, sans-serif; font-size: 11px; margin: 0; padding: 10px; width: 100%%; }\n')];
        report_string = [ report_string, sprintf('.header { margin-bottom: 8px; }\n')];
        report_string = [ report_string, sprintf('h3 { margin: 0 0 3px 0; font-size: 14px; }\n')];
        report_string = [ report_string, sprintf('h4 { margin: 0; font-size: 12px; color: #333; }\n')];
        
        if isIncludeImage
            report_string = [ report_string, sprintf('.image-container { display: flex; justify-content: center; align-items: center; gap: 30px; margin: 8px 0 12px 0; height: 285px; }\n')];
            report_string = [ report_string, sprintf('.image-card { height: 285px; width: auto; max-width: 45%%; object-fit: contain; border: none; background: transparent; }\n')];
        end
        
        report_string = [ report_string, sprintf('.table-wrapper { width: 100%%; padding: 1px 2px; overflow: visible; box-sizing: border-box; }\n')];
        report_string = [ report_string, sprintf('table { border-collapse: collapse; width: calc(100%% - 4px); margin: 5px auto 0 auto; font-size: 10px; table-layout: auto; }\n')];
        report_string = [ report_string, sprintf('th, td { border: 1.2px solid #111; padding: 6px 8px; text-align: center; }\n')];
        report_string = [ report_string, sprintf('th:first-child { text-align: left; padding-left: 10px; }\n')];
        report_string = [ report_string, sprintf('th { background-color: #f2f2f2; font-weight: bold; }\n')];
        report_string = [ report_string, sprintf('.active { background-color: #90EE90 !important; font-weight: bold; -webkit-print-color-adjust: exact; print-color-adjust: exact; }\n')];
        report_string = [ report_string, sprintf('</style>\n</head>\n<body>\n')];
        
        report_string = [ report_string, sprintf('<div class="header">\n')];
        report_string = [ report_string, sprintf('<h3>%s - %s</h3>\n', sInputs(iFile).SubjectName, title_montage)];
        report_string = [ report_string, sprintf('</div>\n')];
        
        if isIncludeImage
            report_string = [ report_string, sprintf('<div class="image-container">\n')];
            if ~isempty(image_path_1)
                img1_b64 = image_to_base64(image_path_1);
                report_string = [ report_string, sprintf('<img class="image-card" src="%s" alt="Image 1">\n', img1_b64)];
            end
            if ~isempty(image_path_2)
                img2_b64 = image_to_base64(image_path_2);
                report_string = [ report_string, sprintf('<img class="image-card" src="%s" alt="Image 2">\n', img2_b64)];
            end
            report_string = [ report_string, sprintf('</div>\n')];
        end
        
        report_string = [ report_string, sprintf('<h3>Source-detector matrix (Get channels info process NIRSTORM) : </h3>\n')];
        report_string = [ report_string, sprintf('<div class="table-wrapper">\n')];
        report_string = [ report_string, sprintf('<table>\n')];
        for row = 1:size(grid_text, 1)
            report_string = [ report_string, sprintf('  <tr>\n')];
            for col = 1:size(grid_text, 2)
                val = grid_text{row, col};
                cls = grid_class{row, col};
                
                if row == 1 || col == 1
                    report_string = [ report_string, sprintf('<th>%s</th>\n', val)];
                else
                    if strcmp(cls, 'active')
                        report_string = [ report_string, sprintf('<td class="active">%s</td>\n', val)];
                    else
                        report_string = [ report_string, sprintf('<td>%s</td>\n', val)];
                    end
                end
            end
            report_string = [ report_string, sprintf('</tr>\n')];
        end
        
        report_string = [ report_string, sprintf('</table>\n')];
        report_string = [ report_string, sprintf('<br /></div>')];
        report_string = [ report_string, sprintf('<h3>Notes: </h3>\n')];
        report_string = [ report_string, sprintf('</body>\n</html>')];
        
        if isSaveToFile
            temp_html = [tempname, '.html'];
            fid_temp = fopen(temp_html, 'w', 'native', 'UTF-8');
            if fid_temp == -1
                bst_error(['Impossible to create temporary file: ', temp_html]);
                return;
            end
            fprintf(fid_temp, '%s', report_string);
            fclose(fid_temp);
            
            a4_width = 1123;
            a4_height = 794;
            
            fig = uifigure('Position', [-5000, -5000, a4_width, a4_height], 'Visible', 'on');
            
            drawnow;
            
            hHtml = uihtml(fig, 'HTMLSource', temp_html, 'Position', [1, 1, a4_width, a4_height]);
            
            drawnow;
            pause(2.0);

            exportapp(fig, pdf_filename);
            
            close(fig);
            if exist(temp_html, 'file')
                delete(temp_html);
            end
            
            if exist(pdf_filename, 'file')
                open(pdf_filename);
            end
        else
            web(sprintf('text://%s', report_string));
        end
    end
end
%% ===== HELPER : IMAGE TO BASE64 DATA URI =====
function b64Uri = image_to_base64(imgPath)
    b64Uri = '';
    if isempty(imgPath) || ~exist(imgPath, 'file')
        return;
    end
    [~, ~, ext] = fileparts(imgPath);
    ext = lower(strrep(ext, '.', ''));
    if strcmp(ext, 'jpg'), ext = 'jpeg'; end
    
    fid = fopen(imgPath, 'rb');
    if fid == -1, return; end
    bytes = fread(fid, inf, 'uint8=>uint8');
    fclose(fid);
    
    encoded = matlab.net.base64encode(bytes);
    b64Uri = sprintf('data:image/%s;base64,%s', ext, encoded);
end

function d = euc_dist(p1, p2)
    d = sqrt(sum((p1 - p2).^2));
end