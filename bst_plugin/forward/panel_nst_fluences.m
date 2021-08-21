function varargout = panel_nst_fluences(varargin)
% panel_nst_fluences Edit ROI for fluences computation. 
% USAGE:  bstPanel = panel_nst_OM('CreatePanel', OPTIONS)           : Call from the interactive interface
%         bstPanel = panel_nst_OM('CreatePanel', sProcess, sFiles)  : Call from the process editor
%                s = panel_nst_OM('GetPanelContents')

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2020 University of Southern California & McGill University
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
% Authors: Francois Tadel, 2020

eval(macro_method);
end


%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(sProcess, sFiles) %#ok<DEFNU>
    panelName = 'FluenceOptions';
    % Java initializations
    import java.awt.*;
    import javax.swing.*;
    import org.brainstorm.list.*;

    ctrl = struct();
    % GUI CALL:  panel_femcond('CreatePanel', OPTIONS)
    if (nargin == 1)
        OPTIONS = sProcess;
    % PROCESS CALL:  panel_femcond('CreatePanel', sProcess, sFiles)
    else
        if isempty(sFiles.ChannelFile)
            sSubject = bst_get('Subject', sProcess.options.subjectname.Value);
            
            OPTIONS.fromMontage = 0;
            if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
                    error('No available Cortex and Head surface for this subject.');
            end    
            OPTIONS.HeadFile = file_fullpath(sSubject.Surface(sSubject.iScalp).FileName);
            OPTIONS.CortexFile = file_fullpath(sSubject.Surface(sSubject.iCortex).FileName);
        else
            OPTIONS.fromMontage = 1;
            OPTIONS.ChannelFile = sFiles.ChannelFile;
            OPTIONS.SubjectName = sFiles.SubjectName;
        end
    end
    
    if OPTIONS.fromMontage 
        ctrl.ChannelFile = sFiles.ChannelFile;
        ctrl.SubjectName = OPTIONS.SubjectName;
    else
        % Load head atlases
        AtlasHead = load(OPTIONS.HeadFile,  'Atlas', 'iAtlas');
        ctrl.HeadAtlasName = {AtlasHead.Atlas.Name};
        % Load cortex atlases
        AtlasCortex = load(OPTIONS.CortexFile, 'Atlas', 'iAtlas');
        ctrl.CortexAtlasName = {AtlasCortex.Atlas.Name};
    end    
    
    % ==== FRAME STRUCTURE ====
    % Create main panel
    jPanelMain = gui_component('Panel');
    % Left panel
    jPanelLeft = gui_river(); 
    jPanelMain.add(jPanelLeft, BorderLayout.CENTER);
    % Right panel
    jPanelRight = gui_river(); 
    jPanelMain.add(jPanelRight, BorderLayout.EAST);
    
    % === PANEL: SEARCH SPACE ===
    jPanelSearchSpace = gui_river([2,2], [3,3,3,3], 'Search space');
        jGroupRadio = ButtonGroup();
        jRadioLayerMontage = gui_component('radio', jPanelSearchSpace, [], 'Montage', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
        jRadioLayerCortex = gui_component('radio', jPanelSearchSpace, [], 'Cortex', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
        jRadioLayerHead = gui_component('radio', jPanelSearchSpace, [], 'Head', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
        if OPTIONS.fromMontage
            jRadioLayerMontage.setSelected(1);
            jRadioLayerCortex.setEnabled(0);
            jRadioLayerHead.setEnabled(0);
        else    
            jRadioLayerMontage.setEnabled(0);
            jRadioLayerCortex.setSelected(1);
        end
        % Register handles
        ctrl.jRadioLayerMontage = jRadioLayerMontage;
        ctrl.jRadioLayerCortex  = jRadioLayerCortex;
        ctrl.jRadioLayerHead    = jRadioLayerHead;
    jPanelLeft.add('hfill', jPanelSearchSpace);
    
    % === PANEL: SCOUTS ===
    % Create panel
    jPanelScouts = gui_river([2,2], [2,7,-10,7], 'Scouts');
        % Create list
        jList = java_create('javax.swing.JList');
        jList.setLayoutOrientation(jList.HORIZONTAL_WRAP);
        jList.setVisibleRowCount(-1);
        % Title
        gui_component('label', jPanelScouts, [], ' Select scouts:', [], [], [], []);
        % Horizontal glue
        gui_component('label', jPanelScouts, 'hfill', ' ', [], [], [], []);
        % Atlas selection box
        jCombo = gui_component('combobox', jPanelScouts, 'right', [], [], [], []);

        % Create scroll panel
        jPanelScouts.add(jList);
        jScroll = javax.swing.JScrollPane(jList);
        jPanelScouts.add('br hfill vfill', jScroll);
        % Set preferred size for the container
        prefPanelSize = java_scaled('dimension', 250,250);
        jPanelScouts.setPreferredSize(prefPanelSize)

        % Extent label
        jExtentTitle = gui_component('label', jPanelScouts, 'br', 'Extent of scalp projection (cm):', [], [], [], []);
        jExtent = gui_component('text', jPanelScouts, 'hfill', '(Default)', [], [], [], []);

        % Register handles
        ctrl.jCombo = jCombo;
        ctrl.jList  = jList ;
        ctrl.jExtent = jExtent;
        ctrl.jPanelScouts = jPanelScouts;
    jPanelLeft.add('br hfill vfill', jPanelScouts);
    
    
    % === PANEL: DESCRIPTION ===
    jPanelRef = gui_river([2,2], [3,5,3,5], 'Description');
        gui_component('label', jPanelRef, '', [...
            '<html>This process uses native MEX version of Monte Carlo eXtreme (MCX)<br />' ...
            'to solve the fluences of each optode. For more details plese refer to:<br /><br />' ...
            'Qianqian Fang and David A. Boas, "Monte Carlo Simulation of Photon<br />'  ...
            'Migrationin 3D Turbid Media Accelerated by Graphics Processing Units".<br />', ...
            'Opt. Express, vol. 17, issue 22, pp. 20178-20190 (2009)</b><br /><br />', ...
            'For technical details please refer to mcx homepage: http://mcx.space'], [],[],[],[]);
    jPanelRight.add('hfill', jPanelRef);

    % === PANEL: OTHER STUFF ====
    jPanelOther = gui_river([2,2], [3,5,3,5], 'Description');
        gui_component('label', jPanelOther, '', 'Any other stuff here', [], [], [], []);
    jPanelRight.add('br hfill', jPanelOther);
    
    % ===== VALIDATION BUTTONS =====
    % Separator
    jPanelSep = java_create('javax.swing.JPanel');
    jPanelSep.add(JLabel(' '));
    jPanelRight.add('br hfill vfill', jPanelSep);
    % Cancel
    jCancelBtn = JButton('Cancel');
    java_setcb(jCancelBtn, 'ActionPerformedCallback', @ButtonCancel_Callback);
    jPanelRight.add('br right', jCancelBtn);
    % Save
    jSearchBtn = JButton('OK');
    java_setcb(jSearchBtn, 'ActionPerformedCallback', @ButtonOk_Callback);
    jPanelRight.add(jSearchBtn);

    
    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Create the BstPanel object that is returned by the function
    bstPanelNew = BstPanel(panelName, jPanelMain, ctrl);    
    % Redraw panel
    UpdatePanel();
    
    
    
    
%% =================================================================================
%  === LOCAL CALLBACKS  ============================================================
%  =================================================================================
    %% ===== BUTTON: CANCEL =====
    function ButtonCancel_Callback(varargin)
        % Close panel
        gui_hide(panelName); % Close panel
        bst_mutex('release', panelName); % Release the MUTEX
    end

    %% ===== BUTTON: OK =====
    function ButtonOk_Callback(varargin)
        % Release mutex and keep the panel opened
        bst_mutex('release', panelName);
    end


    %% ===== UPDATE PANEL ======
    function UpdatePanel()
        if jRadioLayerCortex.isSelected()
            jPanelScouts.setVisible(1);
            jExtentTitle.setVisible(1);
            jExtent.setVisible(1);
            UpdateScoutList(AtlasCortex.Atlas, AtlasCortex.iAtlas);
        elseif jRadioLayerHead.isSelected() 
            jExtentTitle.setVisible(0);
            jExtent.setVisible(0);
            jPanelScouts.setVisible(1);
            UpdateScoutList(AtlasHead.Atlas, AtlasHead.iAtlas);
        elseif jRadioLayerMontage.isSelected()
            jExtentTitle.setVisible(0);
            jExtent.setVisible(0);
            jPanelScouts.setVisible(0);
        end
        
%         % Get panel
%         [bstPanel iPanel] = bst_get('Panel', 'FluenceOptions');
%         container = get(bstPanel, 'container');
%         % Re-pack frame
%         if ~isempty(container)
%             jFrame = container.handle{1};
%             if ~isempty(jFrame)
%                 jFrame.pack();
%             end
%         end
    end


    %% ===== UPDATE SCOUT LIST =====
    function UpdateScoutList(Atlas, iAtlas)
        % Empty the atlas list
        jCombo.removeAllItems();
        % List atlases
        AtlasList = cell(length(Atlas),2);
        for i = 1:length(Atlas)
            AtlasList{i,1} = Atlas(i).Name;
            if ~isempty(Atlas(i).Scouts)
                AtlasList{i,2} = {Atlas(i).Scouts.Label};
            else
                AtlasList{i,2} = [];
            end
            % Add item to the atlas list
            jCombo.addItem(BstListItem(Atlas(i).Name, '', Atlas(i).Name, []));
        end
        
        % Selected atlas
        SelAtlasName = Atlas(iAtlas).Name;
        if ~isempty(SelAtlasName)
            iAtlasList = find(strcmpi({Atlas.Name}, SelAtlasName));
            if isempty(iAtlasList)
                iAtlasList = 1;
            elseif (length(iAtlasList) > 1)
                disp('BST> Error: Two atlases have the same name, you should rename one of them.');
                iAtlasList = iAtlasList(1);
            end
        else
            iAtlasList = 1;
        end
        % Select default atlas
        jCombo.setSelectedIndex(iAtlasList - 1);

        % Set current atlas
        AtlasSelection_Callback(AtlasList, jCombo, jList, []);
        % Set callbacks
        java_setcb(jCombo, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasList, jCombo, jList, ev));
    end
end


%% =================================================================================
%  === EXTERNAL CALLBACKS  =========================================================
%  =================================================================================
%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    % Get panel controls handles
    ctrl = bst_get('PanelControls', 'FluenceOptions');
    if isempty(ctrl)
        s = [];
        return; 
    end
    if ctrl.jRadioLayerCortex.isSelected()
        s.surface = 'cortex';
        s.ROI     = strtrim(char(ctrl.jList.getSelectedValue.getName()));
        s.Atlas   = ctrl.CortexAtlasName(ctrl.jCombo.getSelectedIndex()+1);
        s.Extent  = str2double(ctrl.jExtent.getText);
    elseif ctrl.jRadioLayerHead.isSelected() 
        s.surface = 'head';
        s.ROI     = strtrim(char(ctrl.jList.getSelectedValue.getName()));
        s.Atlas   = ctrl.HeadAtlasName(ctrl.jCombo.getSelectedIndex()+1);
    elseif ctrl.jRadioLayerMontage.isSelected()
        s.surface = 'montage';
        s.ChannelFile = ctrl.ChannelFile;
        s.SubjectName =  ctrl.SubjectName;
    end
end


%% ===== OPTIONS: ATLAS SELECTION CALLBACK =====
function AtlasSelection_Callback(AtlasList, jCombo, jList, ev)
    import org.brainstorm.list.*;
    % Skip deselected event
    if ~isempty(ev) && (ev.getStateChange() ~= ev.SELECTED)
        return;
    end
    % Get current atlas
    iAtlasList = jCombo.getSelectedIndex() + 1;
    if (iAtlasList <= 0)
        return;
    end
    % Get current scouts
    ScoutNames = AtlasList{iAtlasList,2};
    % Temporality disables JList selection callback
    jListCallback_bak = java_getcb(jList, 'ValueChangedCallback');
    java_setcb(jList, 'ValueChangedCallback', []);
    % Create a list of the existing scouts
    listModel = java_create('javax.swing.DefaultListModel');
    for iScout = 1:length(ScoutNames)
        listModel.addElement(BstListItem(ScoutNames{iScout}, '', [' ' ScoutNames{iScout} ' '], iScout));
    end
    jList.setModel(listModel);

    % If there are scouts in this model
    if ~isempty(ScoutNames)
        iSelScouts = 1:length(ScoutNames);
        % Select scouts in the list
        jList.setSelectedIndices(iSelScouts - 1);
        % Save the current selection of scouts (to have the correct list of scouts if the user does not change the selection)
        if jList.isEnabled()
            newList = {AtlasList{iAtlasList,1}, ScoutNames(iSelScouts)};
        else
            newList = {};
        end
        %SetOptionValue(iProcess, optName, newList);
    end
    % Restore JList callback
    java_setcb(jList, 'ValueChangedCallback', jListCallback_bak);
end



