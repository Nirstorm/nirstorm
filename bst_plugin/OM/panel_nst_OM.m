function varargout = panel_nst_OM(varargin)
% panel_nst_OM Panel for the optimal montage
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
% Authors: Francois Tadel, 2020, Edouard Delaire, 2025

eval(macro_method);
end


%% ===== CREATE PANEL =====
function [bstPanelNew, panelName] = CreatePanel(sProcess, sFiles) 
    panelName = 'OMOptions';
    % Java initializations
    import java.awt.*;
    import javax.swing.*;
    import org.brainstorm.list.*;


    ctrl = struct();
    % GUI CALL:  panel_nst_OM('CreatePanel', OPTIONS)
    if (nargin == 1)
        OPTIONS = sProcess;
    % PROCESS CALL:  panel_nst_OM('CreatePanel', sProcess, sFiles)
    else
        sSubject = bst_get('Subject', sProcess.options.subjectname.Value);
        OPTIONS.SubjectName = sProcess.options.subjectname.Value;
        if isempty(sSubject.iCortex) || isempty(sSubject.iScalp)
                error('No available Cortex and Head surface for this subject.');
        end    
        OPTIONS.HeadFile = file_fullpath(sSubject.Surface(sSubject.iScalp).FileName);
        OPTIONS.CortexFile = file_fullpath(sSubject.Surface(sSubject.iCortex).FileName);
    end
    
    OPTIONS = struct_copy_fields(OPTIONS,  getDefaultOptions(), 0);

    if isfield(sProcess.options, 'fluencesCond') && isfield(sProcess.options.fluencesCond,'Value') && ~isempty(sProcess.options.fluencesCond.Value)
        OPTIONS = struct_copy_fields(OPTIONS,  sProcess.options.fluencesCond.Value, 1);
    end

    ctrl.SubjectName = OPTIONS.SubjectName;
    % Load head atlases
    AtlasHead = load(OPTIONS.HeadFile,  'Atlas', 'iAtlas');
    ctrl.HeadAtlasName = {AtlasHead.Atlas.Name};
    
    if ~isempty(OPTIONS.Atlas_head)
        AtlasHead.iAtlas = find(strcmp({AtlasHead.Atlas.Name}, OPTIONS.Atlas_head));
        AtlasHead.iScout = find(strcmp({AtlasHead.Atlas(AtlasHead.iAtlas).Scouts.Label}, OPTIONS.ROI_head));
    else
        AtlasHead.iScout = 1;
    end

    % Load cortex atlases
    AtlasCortex = load(OPTIONS.CortexFile, 'Atlas', 'iAtlas');
    ctrl.CortexAtlasName = {AtlasCortex.Atlas.Name};
    if ~isempty(OPTIONS.Atlas_cortex)
        AtlasCortex.iAtlas = find(strcmp({AtlasCortex.Atlas.Name}, OPTIONS.Atlas_cortex));
        AtlasCortex.iScout = find(strcmp({AtlasCortex.Atlas(AtlasCortex.iAtlas).Scouts.Label}, OPTIONS.ROI_cortex));
    else
        AtlasCortex.iScout = 1;
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
    % === PANEL: SCOUTS ===
    % Create panel
    jPanelScoutsCortex = gui_river([2,2], [2,7,-10,7], 'Cortical scout (target ROI)');
    
    
    % Create list
    jListCortex = java_create('javax.swing.JList');
    jListCortex.setLayoutOrientation(jListCortex.HORIZONTAL_WRAP);
    jListCortex.setVisibleRowCount(-1);
    
    selectionModel = jListCortex.getSelectionModel();
    set(selectionModel, 'ValueChangedCallback', @(e,v)checkSelectionROI());
    
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
    
    

        
    % === PANEL: SCOUTS ===
    % Create panel
    jPanelUseDefault = gui_river([2,2], [3,3,3,3], '');
    jUseDefaultSpace = gui_component('checkbox', jPanelUseDefault, 'br', 'Use default search space', [], [], @(h,ev)UpdatePanel(), []);
    jUseDefaultSpace.setSelected(isempty(OPTIONS.ROI_head));    
    gui_component('label', jPanelUseDefault, 'br', '', [], [], [], []);
    jPanelLeft.add('br hfill vfill', jPanelUseDefault);
    ctrl.jUseDefaultSpace = jUseDefaultSpace;
    
    
    % Extent label
    jExtentTitle = gui_component('label', jPanelUseDefault, 'br', 'Extent of scalp projection:', [], [], [], []);
    jExtent = gui_component('text', jPanelUseDefault, 'hfill', num2str(OPTIONS.Extent), [], [], @(h,ev)checkScalpProjection(), []);
    jExtentTitle2 = gui_component('label', jPanelUseDefault, 'hfill', 'cm', [], [], [], []);
    ctrl.jExtent= jExtent;
        
        
    jPanelScoutsHead = gui_river([2,2], [2,7,-10,7], 'Head scout (search space)');
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
    
    % === PANEL: Objective function ====
    jPanelMontage = gui_river([2,2], [3,5,3,5], 'Objective function');
    
    gui_component('label', jPanelMontage, 'br', 'Criteria to optimize:  ', [], [], [], []);
    jBtnSensitivity = gui_component('checkbox', jPanelMontage, '', 'Sensitivity', [], [], @(h,ev)UpdatePanel(), []);
    jBtnSensitivity.setSelected(1);
    jBtnSensitivity.setEnabled(0);    
 
    gui_component('label', jPanelMontage, '', '   ', [], [], [], []);
    jPanelLeft.add('br hfill vfill', jPanelMontage);
    
    jBtnCoverage = gui_component('checkbox', jPanelMontage, '', 'Coverage', [], [], @(h,ev)UpdatePanel(), []); 
    gui_component('label', jPanelMontage, '', '  ', [], [], [], []);
    jPanelLeft.add('br hfill vfill', jPanelMontage); 
    jBtnCoverage.setSelected(OPTIONS.include_coverage); 
    jBtnCoverage.setEnabled(1);
    
    ctrl.jBtnSensitivity  = jBtnSensitivity;
    ctrl.jBtnCoverage     = jBtnCoverage;
    
    jPanelOptionCoverage = gui_river([2,2], [3,5,3,5], '');

    gui_component('label', jPanelOptionCoverage, 'br', 'Lambda values (Range):  Min = ', [], [], [], []);
    jLambda_min = gui_component('text', jPanelOptionCoverage, 'hfill', num2str(OPTIONS.lambda_coverage(1)), [], [], @(h,ev)checkLimitLambda(), []);
    gui_component('label', jPanelOptionCoverage, '', ' ; Step  = ', [], [], [], []);
    jLambda_step = gui_component('text', jPanelOptionCoverage, 'hfill', num2str(OPTIONS.lambda_coverage(2)), [], [], @(h,ev)checkLimitLambda(), []);
    gui_component('label', jPanelOptionCoverage, '', ' ; Max  = ', [], [], [], []);
    jLambda_max = gui_component('text', jPanelOptionCoverage, 'hfill', num2str(OPTIONS.lambda_coverage(3)), [], [], @(h,ev)checkLimitLambda(), []);
    
    gui_component('label', jPanelOptionCoverage, 'hfill', ' ', [], [], [], []);
    
    jPanelMontage.add('br hfill', jPanelOptionCoverage);

    ctrl.jLambda_min  = jLambda_min;
    ctrl.jLambda_step = jLambda_step;
    ctrl.jLambda_max  = jLambda_max;
    
    gui_component('label', jPanelMontage, 'br', 'Folder for weight table:', [], [], [], []);
    jWeightFolder = gui_component('text', jPanelMontage, 'hfill', OPTIONS.outputdir, [], [], @(h, ev)checkFolder(), []);
    ctrl.jWeightFolder   = jWeightFolder;

    jPanelRight.add('br hfill', jPanelMontage);
    
    % === PANEL: Montage information ====
    jPanelMontage = gui_river([2,2], [3,5,3,5], 'Montage');
    gui_component('label', jPanelMontage, 'br', 'Number of sources:', [], [], [], []);
    jSources = gui_component('text', jPanelMontage, 'hfill', num2str(OPTIONS.nb_sources), [], [], [], []);
    ctrl.jSources = jSources;
    java_setcb(jSources, 'KeyTypedCallback', @(h,ev)isPositive(jSources, 'sources'));

    
    gui_component('label', jPanelMontage, 'br', 'Number of detectors:', [], [], [], []);
    jDetectors = gui_component('text', jPanelMontage, 'hfill', num2str(OPTIONS.nb_detectors), [], [], [], []);
    ctrl.jDetectors = jDetectors;
    java_setcb(jDetectors, 'KeyTypedCallback', @(h,ev)isPositive(jDetectors, 'detectors'));
    
    
    gui_component('label', jPanelMontage, 'br', 'Number of adjacence:', [], [], [], []);
    jAdjacent = gui_component('text', jPanelMontage, 'hfill', num2str(OPTIONS.nAdjacentDet), [], [], [], []);
    ctrl.jAdjacent = jAdjacent;
    java_setcb(jAdjacent, 'KeyTypedCallback', @(h,ev)isPositive(jAdjacent, 'adjacence'));
    
        
    gui_component('label', jPanelMontage, 'br', 'Range of optodes distance', [], [], [], []);
    jSepOptodeMin = gui_component('text', jPanelMontage, 'hfill', num2str(OPTIONS.sep_optode(1)), [], [], @(h, ev)checkLimitSepOptode(), []);       
    gui_component('label', jPanelMontage, '', ' - ', [], [], [], []);
    jSepOptodeMax = gui_component('text', jPanelMontage, 'hfill', num2str(OPTIONS.sep_optode(2)), [], [], @(h, ev)checkLimitSepOptode(), []);
    gui_component('label', jPanelMontage, 'hfill', ' mm', [], [], [], []);
    ctrl.jSepOptodeMin = jSepOptodeMin;
    ctrl.jSepOptodeMax = jSepOptodeMax;

    gui_component('label', jPanelMontage, 'br', 'Minimum source detector distance:', [], [], [], []);
    jSepmin_SD = gui_component('text', jPanelMontage, 'hfill', num2str(OPTIONS.sepmin_SD), [], [], @(h, ev)checkLimitSepOptode(), []);    
    gui_component('label', jPanelMontage, 'hfill', ' mm', [], [], [], []);
    ctrl.jSepmin_SD = jSepmin_SD;

    jPanelRight.add('br hfill', jPanelMontage);
    
    % === PANEL: Fluence  ====
    jPanelFluence = gui_river([2,2], [3,5,3,5], 'Fluence information');
    gui_component('label', jPanelFluence, 'br', 'Fluence data source (URL or path):', [], [], [], []);
    jFluenceSource = gui_component('text', jPanelFluence, 'hfill', OPTIONS.data_source, [], [], @(h, ev)checkFluences(), []);
    ctrl.jFluenceSource = jFluenceSource;

    gui_component('label', jPanelFluence, 'br', 'Wavelength (nm)', [], [], [], []);
    jWavelengths = gui_component('text', jPanelFluence, 'hfill', num2str(OPTIONS.wavelengths), [], [], @(h,ev)checkWavelength(), []);
    ctrl.jWavelengths = jWavelengths;

    prefPanelSize = java_scaled('dimension', 800,100);
    jPanelFluence.setPreferredSize(prefPanelSize)
    jPanelRight.add('br hfill', jPanelFluence);
    
    % === PANEL: Output  ====
    jPanelOutput = gui_river([2,2], [3,5,3,5], 'Output');
    gui_component('label', jPanelOutput, 'br', 'Output condition name:', [], [], [], []);
    jOutputCondition = gui_component('text', jPanelOutput, 'hfill', OPTIONS.condition_name, [], [], @(h, ev)condition_name(), []);
    ctrl.jOutputCondition   = jOutputCondition;
    jPanelRight.add('br hfill', jPanelOutput);

     
    % === PANEL: Error  ====
    jPanelError = gui_river([2,2], [3,5,3,5], 'Error message');
    jLabelError = gui_component('label', jPanelError, 'br' , '', [], [], [], []);
    jPanelRight.add('br hfill', jPanelError);

    
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

    % Error control
    errorList = containers.Map();
    
    % ===== PANEL CREATION =====
    % Return a mutex to wait for panel close
    bst_mutex('create', panelName);
    % Create the BstPanel object that is returned by the function
    bstPanelNew = BstPanel(panelName, jPanelMain, ctrl); 


    
    % Redraw panel
    UpdateScoutList(jComboHead, jListHead, AtlasHead.Atlas, AtlasHead.iAtlas, AtlasHead.iScout);
    UpdateScoutList(jComboCortex,jListCortex, AtlasCortex.Atlas, AtlasCortex.iAtlas, AtlasCortex.iScout);
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
        
        if jUseDefaultSpace.isSelected()
            jPanelScoutsHead.setVisible(0);
            jExtentTitle.setVisible(1);
            jExtent.setVisible(1);
            jExtentTitle2.setVisible(1);
        else
            jPanelScoutsHead.setVisible(1);
            jExtentTitle.setVisible(0);
            jExtent.setVisible(0);
            jExtentTitle2.setVisible(0);
        end    
        
        
        if jBtnCoverage.isSelected()
            jPanelOptionCoverage.setVisible(1);
        else
            jPanelOptionCoverage.setVisible(0);
        end
        
        validateButtonOk();
    end

    %% ===== CHECK FIELDS ======
    
    function checkSelectionROI()
        selected_value = jListCortex.getSelectedValues();
        
        if isempty(selected_value)
            errorList('roi_cortex') = 'You must select exactly one region of interest';
        elseif length(selected_value) > 1
            errorList('roi_cortex') = 'You must select exactly one region of interest : Consider merging the ROIs before calling optimal montage. ';
        else
            if isKey(errorList, 'roi_cortex')
                remove(errorList, 'roi_cortex');
            end
        end
        
        validateButtonOk();
    end
    
    function checkLimitLambda()
        min_value = str2double(char(jLambda_min.getText()));
        max_value = str2double(char(jLambda_max.getText()));
        step = str2double(char(jLambda_step.getText()));

        values = [min_value, max_value, step];
        if any(isnan(values)) || any(values < 0) || isempty(min_value:step:max_value)
            errorList('lambda') = 'Check the definition of lambda. ';
        else
            if isKey(errorList, 'lambda')
                remove(errorList, 'lambda');
            end
        end
        validateButtonOk();
    end

    function checkFolder()
        outputdir = char(jWeightFolder.getText());
        if ~isempty(outputdir) && ~isfolder(outputdir)
            errorList('outputdir') = 'Output directory does not exist. ';
        else
            if isKey(errorList, 'outputdir')
                remove(errorList,'outputdir');
            end
        end
        validateButtonOk();
    end

    function isPositive(jFields, name)
        val = str2double(char(jFields.getText()));
        if isnan(val) || val <= 0 || round(val) ~= val
            errorList(name) = ['Number of ' name ' must be a positive integer.'];
        else
            if isKey(errorList, name)
                remove(errorList, name);
            end
        end
        validateButtonOk();
    end

    function checkLimitSepOptode()
        min_value = str2double(char(jSepOptodeMin.getText()));
        max_value = str2double(char(jSepOptodeMax.getText()));
        sepmin = str2double(char(jSepmin_SD.getText()));
        
        values = [min_value, max_value];
        if any(isnan(values)) || any(values < 0) || any(round(values) ~= values) || isempty(min_value:max_value)
            errorList('SepOptode') = 'The optode distance must be be a valid interval (Ex : 15 - 40). ';
        else
            if isKey(errorList, 'SepOptode')
                remove(errorList, 'SepOptode');
            end
        end
        
        if isnan(sepmin) || sepmin > max_value || sepmin < min_value || round(sepmin) ~= sepmin
            errorList('SepMinOpt') = 'Minimum source detector distance must be an integer within the allowed range. ';
        else
            if isKey(errorList, 'SepMinOpt')
                remove(errorList, 'SepMinOpt');
            end
        end
        validateButtonOk();
    end

    function checkWavelength()
        val = str2double(char(jWavelengths.getText()));
        if isnan(val) || val <= 0 || round(val) ~= val
            errorList('Wavelength') = 'Wavelength must be a positive integer.';
        else
            if isKey(errorList, 'Wavelength')
                remove(errorList, 'Wavelength');
            end
        end
        validateButtonOk();
    end

    function checkScalpProjection()
        val = str2double(char(jExtent.getText()));
        if isnan(val) || val <= 0 || round(val) ~= val
            errorList('ScalpProjection') = 'The value for the extent of scalp projection must be a positive integer.';
        else
            if isKey(errorList, 'ScalpProjection')
                remove(errorList, 'ScalpProjection');
            end
        end
        validateButtonOk();
    end

    function checkFluences()
        %TODO : Check if the fluences directory/website selected has fluences in it
        val = char(jFluenceSource.getText());
        if isempty(val)
            errorList('FluenceSource') = 'Fluences data source must not be empty.';
        else
            if isKey(errorList, 'FluenceSource')
                remove(errorList, 'FluenceSource');
            end
        end
        validateButtonOk();
    end

    function condition_name()
        val = char(jOutputCondition.getText());
        if isempty(val)
            errorList('OutputCondition') = 'You must specifiy an output condition name.';
        else
            if isKey(errorList, 'OutputCondition')
                remove(errorList, 'OutputCondition');
            end
        end
        validateButtonOk();
    end

    function validateButtonOk()
        if isempty(errorList)
            jSearchBtn.setEnabled(1);
            jPanelError.setVisible(0);
            jLabelError.setText('ok')

        else
            jSearchBtn.setEnabled(0);  
            jPanelError.setVisible(1);

            error_keys = errorList.keys();
            error_msg = '<html>';
            
            for iError  = 1:length(error_keys)
                error_msg = [error_msg, '- ' errorList(error_keys{iError}), '<br />'];
            end
            error_msg = [error_msg, '</html>'];

            jLabelError.setText(error_msg)
           
        end
    end

    %% ===== UPDATE SCOUT LIST =====
    function UpdateScoutList(jCombo, jList, Atlas, iAtlas, iScout)
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
        
        if ~isempty(iScout)
            jList.setSelectedIndex(iScout - 1);
        end

        % Set callbacks
        java_setcb(jCombo, 'ItemStateChangedCallback', @(h,ev)AtlasSelection_Callback(AtlasList, jCombo, jList, ev));
    end
end


%% =================================================================================
%  === EXTERNAL CALLBACKS  =========================================================
%  =================================================================================

function options = getDefaultOptions()
    options = struct();
    
    %ROI
    options.surface         = 'cortex';
    options.ROI_cortex      = [];
    options.Atlas_cortex    = [];
    
    options.ROI_head        = [];
    options.Atlas_head      = [];

    options.Extent          = 5;
    
    % Objective function
    options.include_coverage = 0;
    options.lambda_coverage  = [0, 1, 1];
    options.outputdir        = '';
    
    % Montage
    options.nb_sources      = 3;
    options.nb_detectors    = 7;
    options.nAdjacentDet    = 7;
    options.sep_optode      = [ 15, 40];
    options.sepmin_SD       =  15;
    
    %Fluences info
    options.wavelengths     = 685;
    options.data_source     = [nst_get_repository_url() '/fluence/'];
    options.exist_weight    = 1;  
    
    %Output
    options.condition_name  = 'planning_optimal_montage';

end


%% ===== GET PANEL CONTENTS =====
function s = GetPanelContents() 
    % Get panel controls handles
    ctrl = bst_get('PanelControls', 'OMOptions');
    if isempty(ctrl)
        s = [];
        return; 
    end
    
    s.surface = 'cortex';
    s.ROI_cortex     = strtrim(char(ctrl.jListCortex.getSelectedValue.getName()));
    s.Atlas_cortex   = ctrl.CortexAtlasName(ctrl.jComboCortex.getSelectedIndex()+1);

    if ctrl.jUseDefaultSpace.isSelected
        s.ROI_head = [];
        s.Atlas_head = [];
        s.Extent = str2double(ctrl.jExtent.getText);
    else    
        s.ROI_head     = strtrim(char(ctrl.jListHead.getSelectedValue.getName()));
        s.Atlas_head   = ctrl.HeadAtlasName(ctrl.jComboHead.getSelectedIndex()+1);
    end
    s.SubjectName =  ctrl.SubjectName;
    
    % Objective function
    s.include_coverage  = ctrl.jBtnCoverage.isSelected();
    s.lambda_coverage   = [str2double(ctrl.jLambda_min.getText), str2double(ctrl.jLambda_step.getText), str2double(ctrl.jLambda_max.getText)];

    s.outputdir = strtrim(char(ctrl.jWeightFolder.getText));
    
    s.nb_sources = str2double(ctrl.jSources.getText);
    s.nb_detectors = str2double(ctrl.jDetectors.getText);
    s.nAdjacentDet = str2double(ctrl.jAdjacent.getText);
    s.sep_optode  = [ str2double(ctrl.jSepOptodeMin.getText), str2double(ctrl.jSepOptodeMax.getText)];
    s.sepmin_SD  = str2double(ctrl.jSepmin_SD.getText);

    s.wavelengths = strtrim(char(ctrl.jWavelengths.getText));
    
    s.condition_name = strtrim(char(ctrl.jOutputCondition.getText));
    s.data_source = strtrim(char(ctrl.jFluenceSource.getText));
    s.exist_weight = 1;  
    s.flag_display = 0;  

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
    if ~isempty(ScoutNames)
        ScoutNames  =  ScoutNames( ~strcmp(ScoutNames, 'FluenceRegion') & ~strcmp(ScoutNames, 'FluenceExclude'));
    end

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



