function varargout = panel_nst_rename_montage(varargin)
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
% Authors: Jean-Eudes Bornert, 2026 ; Edouard Delaire 2026

    eval(macro_method);
end

%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(sProcess, sFiles) 
    panelName = 'RenameMontageOptions';
    import java.awt.*;
    import javax.swing.*;
    
    ctrl = struct();
    
    jPanelMain = gui_component('Panel');
    jPanelRiver = gui_river([5,5], [5,5,5,5], 'New configuration');
    jPanelMain.add(jPanelRiver, BorderLayout.CENTER);
    
    if isempty(sFiles)
        bstPanelNew = BstPanel(panelName, jPanelMain, ctrl);
        return;
    end
    
    sChannel = in_bst_channel(sFiles(1).ChannelFile);
    sources = []; detectors = [];
    for i = 1:length(sChannel.Channel)
        [isrc, idet, ~, ~] = nst_unformat_channel(sChannel.Channel(i).Name);
        sources(end+1) = isrc; 
        detectors(end+1) = idet;
    end
    sources = unique(sources(sources > 0));
    detectors = unique(detectors(detectors > 0));
    
    prev_transf = struct();
    if isfield(sProcess.options, 'montage_panel') && ...
       isfield(sProcess.options.montage_panel, 'Value') && ...
       ~isempty(sProcess.options.montage_panel.Value) && ...
       isfield(sProcess.options.montage_panel.Value, 'transformations')
   
        prev_transf = sProcess.options.montage_panel.Value.transformations;
    end
    
    ctrl.txtSources = struct();
    if ~isempty(sources)
        gui_component('label', jPanelRiver, 'br', '<HTML><B>Sources</B></HTML>', [], [], [], []);
        for i = 1:length(sources)
            src = sources(i);
            lbl = sprintf('S %d  ->  S ', src);
            gui_component('label', jPanelRiver, 'br', lbl, [], [], [], []);
            
            src_key = sprintf('S%d', src);
            
            if isfield(prev_transf, src_key)
                default_val = strrep(prev_transf.(src_key), 'S', '');
            else
                default_val = num2str(src);
            end
            
            jTxt = gui_component('text', jPanelRiver, 'tab', default_val, [], [], [], []);
            ctrl.txtSources.(src_key) = jTxt;
        end
    end
    
    ctrl.txtDetectors = struct();
    if ~isempty(detectors)
        gui_component('label', jPanelRiver, 'br', '<HTML><B>Detectors</B></HTML>', [], [], [], []);
        for i = 1:length(detectors)
            det = detectors(i);
            lbl = sprintf('D %d  ->  D ', det);
            gui_component('label', jPanelRiver, 'br', lbl, [], [], [], []);
            
            det_key = sprintf('D%d', det);
            
            if isfield(prev_transf, det_key)
                default_val = strrep(prev_transf.(det_key), 'D', '');
            else
                default_val = num2str(det);
            end
            
            jTxt = gui_component('text', jPanelRiver, 'tab', default_val, [], [], [], []);
            ctrl.txtDetectors.(det_key) = jTxt;
        end
    end
    
    jPanelButtons = java_create('javax.swing.JPanel');
    jPanelButtons.setLayout(java_create('java.awt.FlowLayout', java.awt.FlowLayout.RIGHT));
    
    jCancelBtn = javax.swing.JButton('Cancel');
    java_setcb(jCancelBtn, 'ActionPerformedCallback', @ButtonCancel_Callback);
    jPanelButtons.add(jCancelBtn);
    
    jOkBtn = javax.swing.JButton('OK');
    java_setcb(jOkBtn, 'ActionPerformedCallback', @ButtonOk_Callback);
    jPanelButtons.add(jOkBtn);
    
    jPanelMain.add(jPanelButtons, BorderLayout.SOUTH);

    bst_mutex('create', panelName);
    
    bstPanelNew = BstPanel(panelName, jPanelMain, ctrl); 

    %% --- Callbacks ---
    function ButtonCancel_Callback(varargin)
        gui_hide(panelName);
        bst_mutex('release', panelName);
    end

    function ButtonOk_Callback(varargin)
        nb_sources = length(sources);
        src_vals = [];
        srcFields = fieldnames(ctrl.txtSources);
        for iSrc = 1:length(srcFields)
            str_val = char(ctrl.txtSources.(srcFields{iSrc}).getText());
            str_val = strtrim(str_val);
            
            % Empty field
            if isempty(str_val)
                java_dialog('error', 'All fields has to be filled. ', 'Error');
                return;
            end
            
            val = str2double(str_val);
            
            % NaN
            if isnan(val)
                java_dialog('error', 'All fields must be numbers only. ', 'Error');
                return;
            end
            
            % Not an integer
            if round(val) ~= val
                java_dialog('error', 'All fields must be integers', 'Error');
                return;
            end
            
            % Out of bounds
            if val < 1 || val > nb_sources
                java_dialog('error', sprintf('All sources number must be between 1 et %d.', nb_sources), 'Error');
                return;
            end
            
            src_vals(end+1) = val;
        end
        
        % Same number
        if length(unique(src_vals)) < length(src_vals)
            java_dialog('error', '2 sources cannot have the same number.', 'Error');
            return;
        end
        
        nb_detectors = length(detectors);
        det_vals = [];
        detFields = fieldnames(ctrl.txtDetectors);
        for iDet = 1:length(detFields)
            str_val = char(ctrl.txtDetectors.(detFields{iDet}).getText());
            str_val = strtrim(str_val);

            % Empty field
            if isempty(str_val)
                java_dialog('error', 'All fields has to be filled. ', 'Error');
                return;
            end
            
            val = str2double(str_val);
            
            % NaN
            if isnan(val)
                java_dialog('error', 'All fields must be numbers only. ', 'Error');
                return;
            end
            
            % Not an integer
            if round(val) ~= val
                java_dialog('error', 'All fields must be integers', 'Error');
                return;
            end

            % Out of bounds
            if val < 1 || val > nb_detectors
                java_dialog('error', sprintf('All detectors number must be between 1 et %d.', nb_detectors), 'Error');
                return;
            end
            
            det_vals(end+1) = val;
        end
        
        % Same number
        if length(unique(det_vals)) < length(det_vals)
            java_dialog('error', '2 detectors cannot have the same number.', 'Error');
            return; 
        end

        bst_mutex('release', panelName);
    end
end

%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() 
    ctrl = bst_get('PanelControls', 'RenameMontageOptions');
    
    if isempty(ctrl)
        s = []; return; 
    end
    
    s = struct();
    s.transformations = struct();
    
    if isfield(ctrl, 'txtSources')
        srcFields = fieldnames(ctrl.txtSources);
        for i = 1:length(srcFields)
            oldName = srcFields{i};
            newVal = char(ctrl.txtSources.(oldName).getText());
            if ~isempty(newVal)
                s.transformations.(oldName) = sprintf('S%s', strtrim(newVal));
            end
        end
    end
    
    if isfield(ctrl, 'txtDetectors')
        detFields = fieldnames(ctrl.txtDetectors);
        for i = 1:length(detFields)
            oldName = detFields{i};
            newVal = char(ctrl.txtDetectors.(oldName).getText());
            if ~isempty(newVal)
                s.transformations.(oldName) = sprintf('D%s', strtrim(newVal));
            end
        end
    end
end