function varargout = panel_nst_search_space(varargin)

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
% Authors: Jean-Eudes Bornert, 2026, Edouard Delaire, 2026

eval(macro_method);
end

%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(sProcess, sFiles) 
    panelName = 'SearchSpaceOptions';
    import java.awt.*;
    import javax.swing.*;
    import org.brainstorm.list.*;

    ctrl = struct();
    if (nargin == 1)
        OPTIONS = sProcess;
    else
        sSubject = bst_get('Subject', sProcess.options.subjectname.Value);
        
        if isempty(sSubject) || isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
                error('No available Cortex and Head surface for this subject.');
        end    
        OPTIONS.HeadFile = file_fullpath(sSubject.Surface(sSubject.iScalp).FileName);
        OPTIONS.CortexFile = file_fullpath(sSubject.Surface(sSubject.iCortex).FileName);
    end
   
    OPTIONS = struct_copy_fields(OPTIONS, getDefaultOptions(), 0);
    if isfield(sProcess.options.fluencesCond,'Value') && ~isempty(sProcess.options.fluencesCond.Value)
        OPTIONS = struct_copy_fields(OPTIONS,  sProcess.options.fluencesCond.Value, 1);
    end

    AtlasHead = load(OPTIONS.HeadFile, 'Atlas');
    sUserScout = AtlasHead.Atlas(find(strcmp({AtlasHead.Atlas.Name},'User scouts' )));
    if ~isempty(sUserScout.Scouts)
        have_fluence_region     = any(strcmp({sUserScout.Scouts.Label},'FluenceRegion'));
        have_fluence_exclude    = any(strcmp({sUserScout.Scouts.Label},'FluenceExclude'));
    else
        have_fluence_region     = 0;
        have_fluence_exclude    = 0;
    end

    AtlasCortex             = load(OPTIONS.CortexFile, 'Atlas', 'iAtlas');
    ctrl.CortexAtlasName    = {AtlasCortex.Atlas.Name};

    if ~isempty(OPTIONS.Atlas)
        AtlasCortex.iAtlas = find(strcmp({AtlasCortex.Atlas.Name}, OPTIONS.Atlas));
        AtlasCortex.iScout = find(strcmp({AtlasCortex.Atlas(AtlasCortex.iAtlas).Scouts.Label}, OPTIONS.ROI));
    else
        AtlasCortex.iScout = 1;
    end

    % ==== FRAME STRUCTURE ====
    jPanelMain = gui_component('Panel');
    
    jPanelScouts = gui_river([2,2], [2,7,7,7], 'Cortical Scouts');
    jList = java_create('javax.swing.JList');
    jList.setLayoutOrientation(jList.HORIZONTAL_WRAP);
    jList.setVisibleRowCount(-1);
    
    gui_component('label', jPanelScouts, [], ' Select scouts:', [], [], [], []);
    gui_component('label', jPanelScouts, 'hfill', ' ', [], [], [], []);
    jCombo = gui_component('combobox', jPanelScouts, 'right', [], [], [], []);

    jPanelScouts.add(jList);
    jScroll = javax.swing.JScrollPane(jList);
    jPanelScouts.add('br hfill vfill', jScroll);
    prefPanelSize = java_scaled('dimension', 250, 250);
    jPanelScouts.setPreferredSize(prefPanelSize);

    jExtentTitle = gui_component('label', jPanelScouts, 'br', 'Extent of scalp projection:', [], [], [], []);
    jExtent = gui_component('text', jPanelScouts, 'hfill', num2str(OPTIONS.Extent), [], [], [], []);
    gui_component('label', jPanelScouts, 'hfill', ' cm', [], [], [], []);
    
    if have_fluence_region && have_fluence_exclude
        gui_component('label', jPanelScouts, 'br', 'Fluence region and Fluence exclude region detected', [], [], [], []);
    else
        if have_fluence_region
            gui_component('label', jPanelScouts, 'br', 'Fluence region detected', [], [], [], []);
        elseif have_fluence_exclude
            gui_component('label', jPanelScouts, 'br', 'Fluence exclude region detected', [], [], [], []);
        end        
    end
    
    ctrl.jCombo = jCombo;
    ctrl.jList  = jList;
    ctrl.jExtent = jExtent;
    ctrl.jPanelScouts = jPanelScouts;
    
    jPanelMain.add(jPanelScouts, BorderLayout.CENTER);

    % ===== VALIDATION BUTTONS =====
    jPanelBottom = gui_river([0,0], [5,5,5,5]);
    
    jCancelBtn = JButton('Cancel');
    java_setcb(jCancelBtn, 'ActionPerformedCallback', @ButtonCancel_Callback);
    jPanelBottom.add('right', jCancelBtn); % Align to the right
    
    jSearchBtn = JButton('OK');
    java_setcb(jSearchBtn, 'ActionPerformedCallback', @ButtonOk_Callback);
    jPanelBottom.add(jSearchBtn);
    
    jPanelMain.add(jPanelBottom, BorderLayout.SOUTH);
    
    % ===== PANEL CREATION =====
    bst_mutex('create', panelName);
    bstPanelNew = BstPanel(panelName, jPanelMain, ctrl); 
    
    UpdateScoutList(AtlasCortex.Atlas, AtlasCortex.iAtlas, AtlasCortex.iScout);
    
    %% ===== LOCAL CALLBACKS =====
    function ButtonCancel_Callback(varargin)
        gui_hide(panelName);
        bst_mutex('release', panelName);
    end

    function ButtonOk_Callback(varargin)
        bst_mutex('release', panelName);
    end

    function UpdateScoutList(Atlas, iAtlas, iScout)
        import org.brainstorm.list.*;
        jCombo.removeAllItems();
        AtlasList = cell(length(Atlas),2);
        for i = 1:length(Atlas)
            AtlasList{i,1} = Atlas(i).Name;
            if ~isempty(Atlas(i).Scouts)
                AtlasList{i,2} = {Atlas(i).Scouts.Label};
            else
                AtlasList{i,2} = [];
            end
            jCombo.addItem(BstListItem(Atlas(i).Name, '', Atlas(i).Name, []));
        end
        
        SelAtlasName = Atlas(iAtlas).Name;
        if ~isempty(SelAtlasName)
            iAtlasList = find(strcmpi({Atlas.Name}, SelAtlasName));
            if isempty(iAtlasList)
                iAtlasList = 1;
            elseif (length(iAtlasList) > 1)
                iAtlasList = iAtlasList(1);
            end
        else
            iAtlasList = 1;
        end
        jCombo.setSelectedIndex(iAtlasList - 1);
        AtlasSelection_Callback(AtlasList, jCombo, jList, []);

        if ~isempty(iScout)
            jList.setSelectedIndex(iScout - 1);
        end
        java_setcb(jCombo, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasList, jCombo, jList, ev));
    end
end

function options = getDefaultOptions()
    options.surface = 'cortex';
    options.Atlas   = [];
    options.ROI     = [];
    options.Extent  = 4;
end

%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() 
    ctrl = bst_get('PanelControls', 'SearchSpaceOptions');
    if isempty(ctrl)
        s = [];
        return; 
    end
    s.surface = 'cortex';
    tmp = ctrl.jList.getSelectedValuesList.toArray();
    s.ROI  = strjoin(arrayfun(@(x) strtrim(char(x.getName())), tmp, 'UniformOutput', 0 ), ',');
    s.Atlas   = ctrl.CortexAtlasName(ctrl.jCombo.getSelectedIndex()+1);
    s.Extent  = str2double(ctrl.jExtent.getText);
end

%% ===== OPTIONS: ATLAS SELECTION CALLBACK =====
function AtlasSelection_Callback(AtlasList, jCombo, jList, ev)
    import org.brainstorm.list.*;
    if ~isempty(ev) && (ev.getStateChange() ~= ev.SELECTED)
        return;
    end
    iAtlasList = jCombo.getSelectedIndex() + 1;
    if (iAtlasList <= 0)
        return;
    end
    ScoutNames = AtlasList{iAtlasList,2};
    if ~isempty(ScoutNames)
        ScoutNames  =  ScoutNames( ~strcmp(ScoutNames, 'FluenceRegion') & ~strcmp(ScoutNames, 'FluenceExclude'));
    end
    jListCallback_bak = java_getcb(jList, 'ValueChangedCallback');
    java_setcb(jList, 'ValueChangedCallback', []);
    
    listModel = java_create('javax.swing.DefaultListModel');
    for iScout = 1:length(ScoutNames)
        listModel.addElement(BstListItem(ScoutNames{iScout}, '', [' ' ScoutNames{iScout} ' '], iScout));
    end
    jList.setModel(listModel);

    if ~isempty(ScoutNames)
        iSelScouts = 1:length(ScoutNames);
        jList.setSelectedIndices(iSelScouts - 1);
    end
    java_setcb(jList, 'ValueChangedCallback', jListCallback_bak);
end