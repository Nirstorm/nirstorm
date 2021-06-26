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
        % ==== GET Atlas INFO ====
        AtlasHead = load(OPTIONS.HeadFile,  'Atlas', 'iAtlas');
        AtlasCortex = load(OPTIONS.CortexFile, 'Atlas', 'iAtlas');
        
        ctrl.CortexAtlasName = {AtlasCortex.Atlas.Name};
        ctrl.HeadAtlasName = {AtlasHead.Atlas.Name};
    end    
    
    % ==== FRAME STRUCTURE ====
    
    jPanelMain = gui_component('Panel');
    jPanelMain.setLayout(GridLayout(1,2)); % 2 Column
    jPanelMain.setMinimumSize(Dimension(1500,100));

    % Left column
    jPanelLeft = java_create('javax.swing.JPanel');
    jPanelLeft.setLayout(java_create('java.awt.GridBagLayout'));
    
    
    c = GridBagConstraints();
    c.fill = GridBagConstraints.VERTICAL;
    c.gridy = 0;
    c.weightx = 0;
    c.weighty = 0;
    c.insets = Insets(3,5,3,5); % NO IDEA OF WHAT THIS IS
    

    jPanelSearchSpace = java_create('javax.swing.JPanel');
    jPanelSearchSpace.setLayout(GridLayout(1,3));
    jGroupRadio = ButtonGroup();
    jRadioLayerMontage = gui_component('radio', [], 'tab', 'Montage', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
    jRadioLayerCortex = gui_component('radio', [], 'tab', 'Cortex', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
    jRadioLayerHead = gui_component('radio', [], 'tab', 'Head', jGroupRadio, [], @(h,ev)UpdatePanel(), []);
    
    jPanelSearchSpace.add(jRadioLayerMontage);
    jPanelSearchSpace.add(jRadioLayerCortex);
    jPanelSearchSpace.add(jRadioLayerHead);
    
    if OPTIONS.fromMontage
        jRadioLayerMontage.setSelected(1);
        jRadioLayerCortex.setEnabled(0);
        jRadioLayerHead.setEnabled(0);
    else    
        jRadioLayerMontage.setEnabled(0);
        jRadioLayerCortex.setSelected(1);
    end
    
    c.gridy = 0;
    jPanelLeft.add(jPanelSearchSpace, c);

   ctrl.jRadioLayerMontage= jRadioLayerMontage;
   ctrl.jRadioLayerCortex  = jRadioLayerCortex;
   ctrl.jRadioLayerHead   = jRadioLayerHead;

   if ~OPTIONS.fromMontage

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

            gui_component('label', jPanelOptCortex, 'br', 'Extent of scalp projection(cm)');
            jExtent = gui_component('text', jPanelOptCortex, 'hfill', '', [], [], [], []);


            % Set callbacks
            java_setcb(jCombo, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasList, jCombo, jList, ev));
            java_setcb(jList,  'ValueChangedCallback', @(h,ev)ScoutSelection_Callback(AtlasList, jCombo, jList, jCheck, ev));
            %if ~isempty(jCheck)
            %    java_setcb(jCheck, 'ActionPerformedCallback', @(h,ev)ScoutSelection_Callback(iProcess, optNames{iOpt}, AtlasList, jCombo, jList, jCheck, []));
            %end
            % Create scroll panel
        
            jPanelOptCortex.add(jList);
            jScroll = javax.swing.JScrollPane(jList);
            jPanelOptCortex.add('br hfill vfill', jScroll);
            % Set preferred size for the container
            prefPanelSize = java_scaled('dimension', 250,250);
            jPanelOptCortex.setPreferredSize(prefPanelSize)
            c.gridy = 1;
            jPanelLeft.add(jPanelOptCortex, c);
   
            ctrl.jCombo =                   jCombo;
            ctrl.jList  =                   jList ;
            ctrl.jExtent=                   jExtent;
            ctrl.jPanelOptCortex=           jPanelOptCortex;
        end
        
        % Add panel to select scout on the head
        AtlasListHead = cell(length(AtlasHead.Atlas),2);
        SelAtlasName = AtlasHead.Atlas(AtlasHead.iAtlas).Name;

        for i = 1:length(AtlasHead.Atlas)
            AtlasListHead{i,1} = AtlasHead.Atlas(i).Name;
            if ~isempty(AtlasHead.Atlas(i).Scouts)
                AtlasListHead{i,2} = {AtlasHead.Atlas(i).Scouts.Label};
            else
                AtlasListHead{i,2} = [];
            end
        end
        % Selected atlas
        if ~isempty(SelAtlasName)
            iAtlasList = find(strcmpi({AtlasHead.Atlas.Name}, SelAtlasName));
            if isempty(iAtlasList)
                iAtlasList = 1;
            elseif (length(iAtlasList) > 1)
                disp('BST> Error: Two atlases have the same name, you should rename one of them.');
                iAtlasList = AtlasListHead(1);
            end
        else
            iAtlasList = 1;
        end
        % If no scouts are available
        jPanelOptHead =  gui_river([6,6], [-5,6,15,6], 'Head Scout');
        if isempty(AtlasListHead)
            gui_component('label', jPanelOptHead, [], '<HTML>No scouts available.');
        else
            % Create list
            jListHead = java_create('javax.swing.JList');
            jListHead.setLayoutOrientation(jListHead.HORIZONTAL_WRAP);
            jListHead.setVisibleRowCount(-1);
            jPanelOptHead.add(jListHead);
            % Confirm selection

            jCheck = [];
            gui_component('label', jPanelOptHead, [], ' Select scouts:', [],[],[],[]);
            isListEnable = 1;

            % Horizontal glue
            gui_component('label', jPanelOptHead, 'hfill', ' ', [],[],[],[]);
            % Atlas selection box
            jComboHead = gui_component('combobox', jPanelOptHead, 'right', [], {AtlasListHead(:,1)}, [], []);
            % Select default atlas
            if ~isempty(iAtlasList) && (iAtlasList >= 1) && (iAtlasList <= size(AtlasListHead,1))
                jComboHead.setSelectedIndex(iAtlasList - 1);
            end
            % Enable/disable controls
            jListHead.setEnabled(isListEnable);
            jComboHead.setEnabled(isListEnable);
            % Set current atlas
            AtlasSelection_Callback(AtlasListHead, jComboHead, jListHead, []);
            drawnow;
            % Set callbacks
            java_setcb(jComboHead, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasListHead, jComboHead, jListHead, ev));
            java_setcb(jListHead,  'ValueChangedCallback', @(h,ev)ScoutSelection_Callback(AtlasListHead, jComboHead, jListHead, jCheck, ev));
            %if ~isempty(jCheck)
            %    java_setcb(jCheck, 'ActionPerformedCallback', @(h,ev)ScoutSelection_Callback(iProcess, optNames{iOpt}, AtlasList, jCombo, jList, jCheck, []));
            %end
            % Create scroll panel
            jScroll = javax.swing.JScrollPane(jListHead);
            jPanelOptHead.add('br hfill vfill', jScroll);
            % Set preferred size for the container
            prefPanelSize = java_scaled('dimension', 250,180);
            jPanelOptHead.setPreferredSize(prefPanelSize)

            c.gridy = 2;
            jPanelLeft.add(jPanelOptHead, c);
            
            ctrl.jComboHead =                jComboHead;
            ctrl.jListHead  =                 jListHead;
            ctrl.jPanelOptHead =             jPanelOptHead;
        end    
   end
    
    % Right column
    jPanelRight = java_create('javax.swing.JPanel');
    jPanelRight.setLayout(java_create('java.awt.GridBagLayout'));
    
    jPanelREF=  gui_river([6,6], [-5,6,15,6],'' );
    gui_component('label', jPanelREF, '', ['<html>This process uses native MEX version of Monte Carlo eXtreme (MCX) to solve the fluences of each optode.<br />' ...
    'For more details plese refer to Qianqian Fang and David A. Boas,<br />'  ...
    '"Monte Carlo Simulation of Photon Migration in 3D Turbid Media Accelerated by Graphics Processing Units".<br />', ...
    'Opt. Express, vol. 17, issue 22, pp. 20178-20190 (2009)</b><br />', ...
    'For technical details please refer to mcx homepage (http://mcx.sourceforge.net/cgi-bin/index.cgi)</html>'], [],[],[],[]);
    
    c.gridy = 1;
    jPanelRight.add(jPanelREF, c);

    jPanelBtn = gui_component('Panel');
    % Search & Cancel buttons
    jPanelBtnRight = java_create('javax.swing.JPanel');
    jSearchBtn = gui_component('Button', [], [], 'OK', [], [], @ButtonOk_Callback);
    jCancelBtn = gui_component('Button', [], [], 'Cancel', [], [], @ButtonCancel_Callback);
    jPanelBtnRight.add(jSearchBtn);
    jPanelBtnRight.add(jCancelBtn);
    jPanelBtn.add(jPanelBtnRight, BorderLayout.EAST);
    
    c.gridy = 2;
    jPanelRight.add(jPanelBtn, c);
    
    
   % ===== PANEL CREATION =====
    
   jPanelMain.add(jPanelLeft);
   jPanelMain.add(jPanelRight);

       
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
        if ~OPTIONS.fromMontage
            if jRadioLayerCortex.isSelected()
                jPanelOptCortex.setVisible(1);
                jPanelOptHead.setVisible(0);
            else
                jPanelOptCortex.setVisible(0);
                jPanelOptHead.setVisible(1);
            end    
        end    
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
    elseif ctrl.jRadioLayerMontage.isSelected()
        s.surface = 'montage';
        s.ChannelFile = ctrl.ChannelFile;
        s.SubjectName =  ctrl.SubjectName;
    elseif ctrl.jRadioLayerHead.isSelected() 
        s.surface = 'head';
        s.ROI     = strtrim(char(ctrl.jListHead.getSelectedValue.getName()));
        s.Atlas   = ctrl.HeadAtlasName(ctrl.jComboHead.getSelectedIndex()+1);
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


