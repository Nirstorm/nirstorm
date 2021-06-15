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
        OPTIONS = sProcess.options.fluencesCond.Value;
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

    AtlasHead = load(OPTIONS.HeadFile, 'Atlas');
    AtlasCortex = load(OPTIONS.CortexFile, 'Atlas');

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
    
    jPanelHeadAtlas =  gui_river([6,6], [-5,6,15,6], 'Head Scout'); 
    % Combobox
    jComboMethodMEG = gui_component('ComboBox', jPanelHeadAtlas, 'tab hfill', [], [], [], @(h,ev)UpdatePanel(), []);
    
    for i=1:length(AtlasHead.Atlas)
        jComboMethodMEG.addItem(BstListItem('', '', AtlasHead.Atlas(i).Name, []));
    end
    jPanelNew.add(jPanelHeadAtlas);  

    jPanelCortexAtlas =  gui_river([6,6], [-5,6,15,6], 'Cortical Scout');
    jCortexAtlas = gui_component('ComboBox', jPanelCortexAtlas, 'tab hfill', [], [], [], [], []);
    for i=1:length(AtlasCortex.Atlas)
        jCortexAtlas.addItem(BstListItem('', '', AtlasCortex.Atlas(i).Name, []));
    end
    jPanelNew.add(jPanelCortexAtlas); 

        % ===== VALIDATION BUTTONS =====
    jPanelValidation = gui_river([10 0], [6 10 0 10]);
        gui_component('Button', jPanelValidation, 'br right', 'Cancel', [], [], @ButtonCancel_Callback, []);
        gui_component('Button', jPanelValidation, [], 'OK', [], [], @ButtonOk_Callback, []);
    jPanelNew.add(jPanelValidation);

    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Create the BstPanel object that is returned by the function
    ctrl = struct('jPanelSearchSpace',         jPanelSearchSpace);
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
        [bstPanel iPanel] = bst_get('Panel', 'FluenceOptions');
        container = get(bstPanel, 'container');
        % Re-pack frame
        if ~isempty(container)
            jFrame = container.handle{1};
            if ~isempty(jFrame)
                jFrame.pack();
            end
        end
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
    % FEM layers
    for i = 1:length(ctrl.jTextCond)
        s.FemCond(i) = str2double(char(ctrl.jTextCond(i).getText()));
        s.isIsotropic(i) = ctrl.jRadioLayerIso(i).isSelected();
    end
    % Anisotropy options
    if ctrl.jRadioMethodEma.isSelected()
        s.AnisoMethod = 'ema';
    elseif ctrl.jRadioMethodEmaVc.isSelected()
        s.AnisoMethod = 'ema+vc';
    elseif ctrl.jRadioMethodSim.isSelected()
        s.AnisoMethod = 'simulated';
    end
    % Simulated: Ratio
    s.SimRatio = str2double(char(ctrl.jTextSimRatio.getText()));
    % Simulated: Constraint method
    if ctrl.jRadioConstrWang.isSelected()
        s.SimConstrMethod = 'wang';
    elseif ctrl.jRadioConstrWolters.isSelected()
        s.SimConstrMethod = 'wolters';
    end
end




