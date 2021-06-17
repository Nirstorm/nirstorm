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

    % GUI CALL:  panel_femcond('CreatePanel', OPTIONS)
    if (nargin == 1)
        OPTIONS = sProcess;
    % PROCESS CALL:  panel_femcond('CreatePanel', sProcess, sFiles)
    else
        % Get FEM files
        sSubject = bst_get('Subject', sProcess.options.subjectname.Value);
        if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
            error('No available Cortex and Head surface for this subject.');
        end    
        OPTIONS.HeadFile = file_fullpath(sSubject.Surface(sSubject.iScalp).FileName);
        OPTIONS.CortexFile = file_fullpath(sSubject.Surface(sSubject.iCortex).FileName);
    end
    
    % ==== GET Atlas INFO ====
    % Load Atlas
    %AtlasHead = load(OPTIONS.HeadFile, 'Atlas');
    % [AtlasCortex, ~, ~] = panel_scout('GetScouts',OPTIONS.CortexFile);

    AtlasHead = load(OPTIONS.HeadFile,  'Atlas', 'iAtlas');
    AtlasCortex = load(OPTIONS.CortexFile, 'Atlas', 'iAtlas');

    % ==== FRAME STRUCTURE ====
    jPanelNew = java_create('javax.swing.JPanel');
    jPanelNew.setLayout(BoxLayout(jPanelNew, BoxLayout.PAGE_AXIS));
    jPanelNew.setBorder(BorderFactory.createEmptyBorder(12,12,12,12));

    
    jPanelSearchSpace =  gui_river([6,6], [-5,6,15,6], 'Compute fluences from'); 
    jRadioSearchSpace = javaArray('javax.swing.JRadioButton', 3);
    jGroupRadio = ButtonGroup();
    jRadioLayerMontage = gui_component('radio', jPanelSearchSpace, 'tab', 'Montage', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
    jRadioLayerCortex = gui_component('radio', jPanelSearchSpace, 'tab', 'Cortex', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
    jRadioLayerHead = gui_component('radio', jPanelSearchSpace, 'tab', 'Head', jGroupRadio, [], @(h,ev)UpdatePanel(), []);

    jRadioLayerHead.setSelected(1);
    jPanelNew.add(jPanelSearchSpace);

    % Add panel to select scout on the cortex 
    AtlasList = cell(length(AtlasCortex.Atlas),2);
    SelAtlasName = AtlasCortex.Atlas(AtlasCortex.iAtlas).Name;
    
    for i = 1:length(AtlasCortex.Atlas)
        AtlasList{i,1} = AtlasCortex.Atlas(i).Name;
        if ~isempty(AtlasCortex.Atlas(i).Scouts)
            AtlasList{i,2} = {AtlasCortex.Atlas(i).Scouts.Label};
        else
            AtlasList{i,2} = [];
        end
    end
    % Selected atlas
    if ~isempty(SelAtlasName)
        iAtlasList = find(strcmpi({AtlasCortex.Atlas.Name}, SelAtlasName));
        if isempty(iAtlasList)
            iAtlasList = 1;
        elseif (length(iAtlasList) > 1)
            disp('BST> Error: Two atlases have the same name, you should rename one of them.');
            iAtlasList = iAtlasList(1);
        end
    else
        iAtlasList = 1;
    end
    jPanelOptCortex =  gui_river([6,6], [-5,6,15,6], 'Cortical Scout');
    % If no scouts are available
    if isempty(AtlasList)
        gui_component('label', jPanelOptCortex, [], '<HTML>No scouts available.');
    else
        % Create list
        jList = java_create('javax.swing.JList');
        jList.setLayoutOrientation(jList.HORIZONTAL_WRAP);
        jList.setVisibleRowCount(-1);
        % Confirm selection

        jCheck = [];
        gui_component('label', jPanelOptCortex, [], ' Select scouts:', [],[],[],[]);
        isListEnable = 1;
        
        % Horizontal glue
        gui_component('label', jPanelOptCortex, 'hfill', ' ', [],[],[],[]);
        % Atlas selection box
        jCombo = gui_component('combobox', jPanelOptCortex, 'right', [], {AtlasList(:,1)}, [], []);
        

        
        % Select default atlas
        if ~isempty(iAtlasList) && (iAtlasList >= 1) && (iAtlasList <= size(AtlasList,1))
            jCombo.setSelectedIndex(iAtlasList - 1);
        end
        
        % Enable/disable controls
        jList.setEnabled(isListEnable);
        jCombo.setEnabled(isListEnable);

        % Set current atlas
        AtlasSelection_Callback(AtlasList, jCombo, jList, []);
        drawnow;
        
        gui_component('label', jPanelOptCortex, 'br', 'Extent of cortical ROI to scalp projection(cm)');
        jExtent = gui_component('text', jPanelOptCortex, 'hfill', '', [], [], [], []);
        
       
        % Set callbacks
        java_setcb(jCombo, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasList, jCombo, jList, ev));
        java_setcb(jList,  'ValueChangedCallback', @(h,ev)ScoutSelection_Callback(AtlasList, jCombo, jList, jCheck, ev));
        %if ~isempty(jCheck)
        %    java_setcb(jCheck, 'ActionPerformedCallback', @(h,ev)ScoutSelection_Callback(iProcess, optNames{iOpt}, AtlasList, jCombo, jList, jCheck, []));
        %end
        % Create scroll panel
        jScroll = javax.swing.JScrollPane(jList);
        jPanelOptCortex.add('br hfill vfill', jScroll);
        % Set preferred size for the container
        prefPanelSize = java_scaled('dimension', 250,180);
        
        jPanelNew.add(jPanelOptCortex); 
    end
    
    
    % Add panel to select scout on the head
    AtlasList = cell(length(AtlasHead.Atlas),2);
    SelAtlasName = AtlasHead.Atlas(AtlasHead.iAtlas).Name;
    
    for i = 1:length(AtlasHead.Atlas)
        AtlasList{i,1} = AtlasHead.Atlas(i).Name;
        if ~isempty(AtlasHead.Atlas(i).Scouts)
            AtlasList{i,2} = {AtlasHead.Atlas(i).Scouts.Label};
        else
            AtlasList{i,2} = [];
        end
    end
    % Selected atlas
    if ~isempty(SelAtlasName)
        iAtlasList = find(strcmpi({AtlasHead.Atlas.Name}, SelAtlasName));
        if isempty(iAtlasList)
            iAtlasList = 1;
        elseif (length(iAtlasList) > 1)
            disp('BST> Error: Two atlases have the same name, you should rename one of them.');
            iAtlasList = iAtlasList(1);
        end
    else
        iAtlasList = 1;
    end
    % If no scouts are available
    jPanelOptHead =  gui_river([6,6], [-5,6,15,6], 'Head Scout');
    if isempty(AtlasList)
        gui_component('label', jPanelOptHead, [], '<HTML>No scouts available.');
    else
        % Create list
        jListHead = java_create('javax.swing.JList');
        jListHead.setLayoutOrientation(jListHead.HORIZONTAL_WRAP);
        jListHead.setVisibleRowCount(-1);
        % Confirm selection

        jCheck = [];
        gui_component('label', jPanelOptHead, [], ' Select scouts:', [],[],[],[]);
        isListEnable = 1;
        
        % Horizontal glue
        gui_component('label', jPanelOptHead, 'hfill', ' ', [],[],[],[]);
        % Atlas selection box
        jComboHead = gui_component('combobox', jPanelOptHead, 'right', [], {AtlasList(:,1)}, [], []);
        % Select default atlas
        if ~isempty(iAtlasList) && (iAtlasList >= 1) && (iAtlasList <= size(AtlasList,1))
            jComboHead.setSelectedIndex(iAtlasList - 1);
        end
        % Enable/disable controls
        jListHead.setEnabled(isListEnable);
        jComboHead.setEnabled(isListEnable);

        % Set current atlas
        AtlasSelection_Callback(AtlasList, jComboHead, jListHead, []);
        drawnow;
        % Set callbacks
        java_setcb(jComboHead, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasList, jComboHead, jListHead, ev));
        java_setcb(jListHead,  'ValueChangedCallback', @(h,ev)ScoutSelection_Callback(AtlasList, jComboHead, jListHead, jCheck, ev));
        %if ~isempty(jCheck)
        %    java_setcb(jCheck, 'ActionPerformedCallback', @(h,ev)ScoutSelection_Callback(iProcess, optNames{iOpt}, AtlasList, jCombo, jList, jCheck, []));
        %end
        % Create scroll panel
        jScroll = javax.swing.JScrollPane(jListHead);
        jPanelOptCortex.add('br hfill vfill', jScroll);
        % Set preferred size for the container
        prefPanelSize = java_scaled('dimension', 250,180);
        
        jPanelNew.add(jPanelOptHead); 
    end    

        % ===== VALIDATION BUTTONS =====
    jPanelValidation = gui_river([10 0], [6 10 0 10]);
        gui_component('Button', jPanelValidation, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
        gui_component('Button', jPanelValidation, [], 'OK', [], [], @ButtonOk_Callback, []);
    jPanelNew.add(jPanelValidation);

    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Create the BstPanel object that is returned by the function
    ctrl = struct( 'jRadioLayerMontage',        jRadioLayerMontage, ...
                    'jRadioLayerCortex',        jRadioLayerCortex , ...
                    'jRadioLayerHead' ,         jRadioLayerHead, ...
                    'jPanelSearchSpace',        jPanelSearchSpace, ...
                    'jCombo',jCombo, ...
                    'jList',jList , ...
                    'jPanelOptCortex',          jPanelOptCortex); %, ...
                  %'jPanelOptHead',          jPanelOptHead);
                  
    ctrl.CortexAtlasName = {AtlasCortex.Atlas.Name};
    
    % Create the BstPanel object that is returned by the function
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);    
    % Redraw panel
    UpdatePanel();
%% =================================================================================
%  === LOCAL CALLBACKS  ============================================================
%  =================================================================================
    %% ===== BUTTON: CANCEL =====
    function ButtonCancel_Callback(varargin)
        % Close panel
        gui_hide(panelName);
    end

    %% ===== BUTTON: OK =====
    function ButtonOk_Callback(varargin)
        % Release mutex and keep the panel opened
        bst_mutex('release', panelName);
    end

    %% ===== UPDATE PANEL ======
    function UpdatePanel()
        
        if jRadioLayerCortex.isSelected()
            jPanelOptCortex.setEnabled(1);
            jCombo.setEnabled(1);
            jList.setEnabled(1);
        else
            jPanelOptCortex.setEnabled(0);
            jCombo.setEnabled(0);
            jList.setEnabled(0);
        end    
%         % FEM Layers
%         isIsotropic = false(1, nLayers);
%         for j = 1:nLayers
%             isIsotropic(j) = jRadioLayerIso(j).isSelected();
%             % jTextCond(j).setEnabled(isIsotropic(j));
%         end
%         % All isotropic: disable anisotropy options
%         isAllIso = all(isIsotropic);
%         jPanelAniso.setVisible(~isAllIso);
%         % Simulated
%         isSim = ~isAllIso && jRadioMethodSim.isSelected();
%         jPanelSim.setVisible(isSim);
%         
        % Get panel
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
end


%% =================================================================================
%  === EXTERNAL CALLBACKS  =========================================================
%  =================================================================================
%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() %#ok<DEFNU>
    % Get panel controls handles
    disp('hello1')
    ctrl = bst_get('PanelControls', 'FluenceOptions');
    if isempty(ctrl)
        s = [];
        return; 
    end
    if ctrl.jRadioLayerCortex.isSelected()
        s.surface = 'cortex';
        s.ROI     = strtrim(char(ctrl.jList.getSelectedValue.getName()));
        s.Atlas   = ctrl.CortexAtlasName(ctrl.jCombo.getSelectedIndex()+1);
        disp('1')
    elseif ctrl.jRadioLayerMontage.isSelected()
        s.surface = 'montage';
        disp('2')
    elseif ctrl.jRadioLayerHead.isSelected() 
        s.surface = 'head';
        disp('13')
    end    
    
    
    
end


%% ===== OPTIONS: ATLAS SELECTION CALLBACK =====
    function AtlasSelection_Callback(AtlasList, jCombo, jList, ev)
        import org.brainstorm.list.*;
        % Skip deselected event
        if ~isempty(ev) && (ev.getStateChange() ~= ev.SELECTED)
            return;
        end
        % Get current process
        %sCurProcess = GlobalData.Processes.Current(iProcess);
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

%% ===== OPTIONS: SCOUT SELECTION CALLBACK =====
function ScoutSelection_Callback(AtlasList, jCombo, jList, jCheck, ev)
    % Cancel temporary selection
    if ~isempty(ev) && ev.getValueIsAdjusting()
        return;
    end
    % Enable/disable jList and jCombo
    if ~isempty(jCheck)
        isChecked = jCheck.isSelected();
    else
        isChecked = 1;
    end
    jList.setEnabled(isChecked);
    jCombo.setEnabled(isChecked);
    % Get current atlas
    iAtlasList = jCombo.getSelectedIndex() + 1;
    if (iAtlasList <= 0)
        return;
    end
%     % If cluster/scout not selected
    if ~isChecked
%         SetOptionValue(iProcess, optName, []);
%     % If not currently editing
    else
        % Get selected clusters
        iSel = jList.getSelectedIndices() + 1;
        % List of new selected scouts
        newList = AtlasList(iAtlasList,:);
        newList{1,2} = AtlasList{iAtlasList,2}(iSel);
        % Set value
        %SetOptionValue(iProcess, optName, newList);
    end
end

%% ===== GET OPTIONS =====
function GridOptions = GetOptions(ctrl)

disp('hello2')

end
