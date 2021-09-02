function varargout = panel_nst_OM(varargin)
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
    panelName = 'OMOptions';
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
        sSubject = bst_get('Subject', sProcess.options.subjectname.Value);
        OPTIONS.SubjectName = sProcess.options.subjectname.Value;
        if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
                error('No available Cortex and Head surface for this subject.');
        end    
        OPTIONS.HeadFile = file_fullpath(sSubject.Surface(sSubject.iScalp).FileName);
        OPTIONS.CortexFile = file_fullpath(sSubject.Surface(sSubject.iCortex).FileName);
    end
    
    ctrl.SubjectName = OPTIONS.SubjectName;
    % Load head atlases
    AtlasHead = load(OPTIONS.HeadFile,  'Atlas', 'iAtlas');
    ctrl.HeadAtlasName = {AtlasHead.Atlas.Name};
    % Load cortex atlases
    AtlasCortex = load(OPTIONS.CortexFile, 'Atlas', 'iAtlas');
    ctrl.CortexAtlasName = {AtlasCortex.Atlas.Name};
        
    
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
    
    % === PANEL: SCOUTS ===
    % Create panel
    jPanelScoutsHead = gui_river([2,2], [2,7,-10,7], 'Scouts');
    % Create list
    jListHead = java_create('javax.swing.JList');
    jListHead.setLayoutOrientation(jListHead.HORIZONTAL_WRAP);
    jListHead.setVisibleRowCount(-1);
    % Title
    gui_component('label', jPanelScoutsHead, [], ' Select scouts:', [], [], [], []);
    % Horizontal glue
    gui_component('label', jPanelScoutsHead, 'hfill', ' ', [], [], [], []);
    % Atlas selection box
    jComboHead = gui_component('combobox', jPanelScoutsHead, 'right', [], [], [], []);
    % Create scroll panel
    jPanelScoutsHead.add(jListHead);
    jScroll = javax.swing.JScrollPane(jListHead);
    jPanelScoutsHead.add('br hfill vfill', jScroll);
    % Set preferred size for the container
    prefPanelSize = java_scaled('dimension', 250,250);
    jPanelScoutsHead.setPreferredSize(prefPanelSize)
    
    % Register handles
    ctrl.jComboHead = jComboHead;
    ctrl.jListHead  = jListHead ;
    ctrl.jPanelScoutsHead = jPanelScoutsHead;
    jPanelLeft.add('br hfill vfill', jPanelScoutsHead);
    

        % === PANEL: SCOUTS ===
    % Create panel
    jPanelScoutsCortex = gui_river([2,2], [2,7,-10,7], 'Scouts');
    % Create list
    jListCortex = java_create('javax.swing.JList');
    jListCortex.setLayoutOrientation(jListCortex.HORIZONTAL_WRAP);
    jListCortex.setVisibleRowCount(-1);
    % Title
    gui_component('label', jPanelScoutsCortex, [], ' Select scouts:', [], [], [], []);
    % Horizontal glue
    gui_component('label', jPanelScoutsCortex, 'hfill', ' ', [], [], [], []);
    % Atlas selection box
    jComboCortex = gui_component('combobox', jPanelScoutsCortex, 'right', [], [], [], []);
    % Create scroll panel
    jPanelScoutsCortex.add(jListCortex);
    jScroll = javax.swing.JScrollPane(jListCortex);
    jPanelScoutsCortex.add('br hfill vfill', jScroll);
    % Set preferred size for the container
    prefPanelSize = java_scaled('dimension', 250,250);
    jPanelScoutsCortex.setPreferredSize(prefPanelSize)
    
    
    % Register handles
    ctrl.jComboCortex = jComboCortex;
    ctrl.jListCortex  = jListCortex ;
    ctrl.jPanelScoutsCortex = jPanelScoutsCortex;
    jPanelLeft.add('br hfill vfill', jPanelScoutsCortex);
    
    
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
    UpdateScoutList(jComboHead,jListHead, AtlasHead.Atlas, AtlasHead.iAtlas);
    UpdateScoutList(jComboCortex,jListCortex, AtlasCortex.Atlas, AtlasCortex.iAtlas);

    
    
    
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
    function UpdateScoutList(jCombo,jList, Atlas, iAtlas)
        import org.brainstorm.list.*;
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
    ctrl = bst_get('PanelControls', 'OMOptions');
    if isempty(ctrl)
        s = [];
        return; 
    end
    
    s.surface = 'cortex';
    s.ROI_cortex     = strtrim(char(ctrl.jListCortex.getSelectedValue.getName()));
    s.Atlas_cortex   = ctrl.CortexAtlasName(ctrl.jComboCortex.getSelectedIndex()+1);

    s.ROI_head     = strtrim(char(ctrl.jListHead.getSelectedValue.getName()));
    s.Atlas_head   = ctrl.HeadAtlasName(ctrl.jComboHead.getSelectedIndex()+1);
    
    s.SubjectName =  ctrl.SubjectName;
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



