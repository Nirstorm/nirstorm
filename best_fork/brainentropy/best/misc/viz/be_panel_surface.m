function varargout = panel_surface(varargin)
% PANEL_SURFACE: Panel to load and plot surfaces.
% 
% USAGE:  bstPanel = panel_surface('CreatePanel')
%                    panel_surface('UpdatePanel')
%                    panel_surface('CurrentFigureChanged_Callback')
%       nbSurfaces = panel_surface('CreateSurfaceList',      jToolbar, hFig)
%                    panel_surface('UpdateSurfaceProperties')
%         iSurface = panel_surface('AddSurface',             hFig, surfaceFile)
%                    panel_surface('RemoveSurface',          hFig, iSurface)
%                    panel_surface('SetSurfaceTransparency', hFig, iSurf, alpha)
%                    panel_surface('SetSurfaceColor',        hFig, iSurf, colorCortex, colorSulci)
%                    panel_surface('ApplyDefaultDisplay')
%           [isOk] = panel_surface('SetSurfaceData',        hFig, iTess, dataType, dataFile, isStat)
%           [isOk] = panel_surface('UpdateSurfaceData',     hFig, iSurfaces)
%                    panel_surface('UpdateOverlayCube',     hFig, iTess)

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
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
% Authors: Francois Tadel, 2008-2012

macro_methodcall;
end


%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel() %#ok<DEFNU>
    panelName = 'Surface';
    % Java initializations
    import java.awt.*;
    import javax.swing.*;
    import org.brainstorm.icon.*;
   
    % Constants
    LABEL_WIDTH    = 30;
    BUTTON_WIDTH   = 50;
    SLIDER_WIDTH   = 20;
    DEFAULT_HEIGHT = 22;
    TB_HEIGHT      = 28;
    % Create panel
    jPanelNew = JPanel(BorderLayout());
    jPanelNew.setBorder([]);

    % ====================================================================
    % ==== TOOLBAR : SURFACES LIST =======================================
    jToolbar = gui_component('Toolbar', jPanelNew, BorderLayout.NORTH);
        % Title
        gui_component('Label', jToolbar, [], '    Surfaces:');
        % Separation panel
        jToolbar.add(JPanel());
        jToolbar.addSeparator(Dimension(10, TB_HEIGHT));
        % Add/Remove button
        TB_SIZE = Dimension(25,25);
        gui_component('ToolbarButton', jToolbar, [], [], {IconLoader.ICON_SURFACE_ADD, TB_SIZE},    'Add a surface',              @ButtonAddSurfaceCallback);
        gui_component('ToolbarButton', jToolbar, [], [], {IconLoader.ICON_SURFACE_REMOVE, TB_SIZE}, 'Remove surface from figure', @ButtonRemoveSurfaceCallback);
        
    % ====================================================================
    % ==== OPTIONS PANEL ==================================================
    % Create OPTIONS scrollpane
    jPanelOptions = gui_river([3,5]);
        % ===== Create panel : surface configuration =====
        jPanelSurfaceOptions = gui_river([1,1], [1,8,1,4], 'Surface options');              
            % Alpha title 
            jLabelAlphaTitle = gui_component('label', jPanelSurfaceOptions, 'br', 'Transp.:');
            % Alpha slider
            jSliderSurfAlpha = JSlider(0, 100, 0);
            jSliderSurfAlpha.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
            java_setcb(jSliderSurfAlpha, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'SurfAlpha'), ...
                                         'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'SurfAlpha'));
            jPanelSurfaceOptions.add('tab hfill', jSliderSurfAlpha);
            % Alpha label
            jLabelSurfAlpha = gui_component('label', jPanelSurfaceOptions, [], '   0%', {JLabel.RIGHT, Dimension(LABEL_WIDTH, DEFAULT_HEIGHT)});
            % Quick preview
            java_setcb(jSliderSurfAlpha, 'StateChangedCallback',  @(h,ev)SliderQuickPreview(jSliderSurfAlpha, jLabelSurfAlpha, 1));

            % Smooth title
            gui_component('label', jPanelSurfaceOptions, 'br', 'Smooth:');
            % Smooth slider 
            jSliderSurfSmoothValue = JSlider(0, 100, 0);
            jSliderSurfSmoothValue.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
            jSliderSurfSmoothValue.setToolTipText('Smooth surface');
            java_setcb(jSliderSurfSmoothValue, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'SurfSmoothValue'), ...
                                               'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'SurfSmoothValue'));
            jPanelSurfaceOptions.add('tab hfill', jSliderSurfSmoothValue);
            % Smooth ALPHA label
            jLabelSurfSmoothValue = gui_component('label', jPanelSurfaceOptions, [], '   0%', {JLabel.RIGHT, Dimension(LABEL_WIDTH, DEFAULT_HEIGHT)});
            % Quick preview
            java_setcb(jSliderSurfSmoothValue, 'StateChangedCallback',  @(h,ev)SliderQuickPreview(jSliderSurfSmoothValue, jLabelSurfSmoothValue, 1));
            
            % Buttons
            jButtonSurfColor = gui_component('button', jPanelSurfaceOptions, 'br center', 'Color', {Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT), Insets(0,0,0,0)}, 'Set surface color', @ButtonSurfColorCallback);
            jButtonSurfSulci = gui_component('toggle', jPanelSurfaceOptions, '',          'Sulci', {Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT), Insets(0,0,0,0)}, 'Show/hide sulci map', @ButtonShowSulciCallback);
            jButtonSurfEdge  = gui_component('toggle', jPanelSurfaceOptions, '',          'Edge',  {Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT), Insets(0,0,0,0)}, 'Show/hide surface triangles', @ButtonShowEdgesCallback);
        jPanelOptions.add('br hfill', jPanelSurfaceOptions);
    
        % ===== Create panel : data description =====
        jPanelDataOptions = gui_river([1,1], [1,8,1,4], 'Data options');
            % Threshold title
            jLabelThreshTitle = gui_component('label', jPanelDataOptions, [], 'Amplitude:');
            % Threshold slider
            jSliderDataThresh = JSlider(0, 100, 50);
            jSliderDataThresh.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
            java_setcb(jSliderDataThresh, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'DataThreshold'), ...
                                          'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'DataThreshold'));
            jPanelDataOptions.add('tab hfill', jSliderDataThresh);
            % Threshold label
            jLabelDataThresh = gui_component('label', jPanelDataOptions, [], '   0%', {JLabel.RIGHT, Dimension(LABEL_WIDTH, DEFAULT_HEIGHT)});
            % Quick preview
            java_setcb(jSliderDataThresh, 'StateChangedCallback',  @(h,ev)SliderQuickPreview(jSliderDataThresh, jLabelDataThresh, 1));
            
            % Min size title
            jLabelSizeTitle = gui_component('label', jPanelDataOptions, 'br', 'Min size:');
            % Min size slider
            sliderSizeVector = [1:9, 10:2:18, 20:5:60, 70:10:100];
            jSliderSize = JSlider(1, length(sliderSizeVector), 1);
            jSliderSize.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
            java_setcb(jSliderSize, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'SizeThreshold'), ...
                                    'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'SizeThreshold'));
            jPanelDataOptions.add('tab hfill', jSliderSize);
            % Min size label
            jLabelSize = gui_component('label', jPanelDataOptions, [], '   1', {JLabel.RIGHT, Dimension(LABEL_WIDTH, DEFAULT_HEIGHT)});
            % Quick preview
            java_setcb(jSliderSize, 'StateChangedCallback',  @(h,ev)SliderQuickPreview(jSliderSize, jLabelSize, 0));
            
            % Alpha title and slider
            jLabelDataAlphaTitle = gui_component('label', jPanelDataOptions, 'br', 'Transp:');           
            jSliderDataAlpha = JSlider(0, 100, 0);
            jSliderDataAlpha.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
            java_setcb(jSliderDataAlpha, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'DataAlpha'), ...
                                         'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'DataAlpha'));
            jPanelDataOptions.add('tab hfill', jSliderDataAlpha);
            % Data alpha label
            jLabelDataAlpha = gui_component('label', jPanelDataOptions, [], '   0%', {JLabel.RIGHT, Dimension(LABEL_WIDTH, DEFAULT_HEIGHT)});
            % Quick preview
            java_setcb(jSliderDataAlpha, 'StateChangedCallback',  @(h,ev)SliderQuickPreview(jSliderDataAlpha, jLabelDataAlpha, 1));
        jPanelOptions.add('br hfill', jPanelDataOptions);
        
        % ===== Create panel : surface resect =====
        jPanelSurfaceResect = gui_river([0,4], [1,8,8,0], 'Resect [X,Y,Z]');
            % === RESECT SLIDERS ===
            % Sub panel
            panelResect = JPanel();
            panelResect.setLayout(BoxLayout(panelResect, BoxLayout.LINE_AXIS));
                % Resect X : Slider 
                jSliderResectX = JSlider(-100, 100, 0);
                jSliderResectX.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
                java_setcb(jSliderResectX, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'ResectX'), ...
                                           'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'ResectX'));
                panelResect.add('hfill', jSliderResectX);   
                % Resect Y : Title and Slider 
                jSliderResectY = JSlider(-100, 100, 0);
                jSliderResectY.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
                java_setcb(jSliderResectY, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'ResectY'), ...
                                           'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'ResectY'));
                panelResect.add('hfill', jSliderResectY);     
                % Resect Z : Title and Slider 
                jSliderResectZ = JSlider(-100, 100, 0);
                jSliderResectZ.setPreferredSize(Dimension(SLIDER_WIDTH, DEFAULT_HEIGHT));
                java_setcb(jSliderResectZ, 'MouseReleasedCallback', @(h,ev)SliderCallback(h, ev, 'ResectZ'), ...
                                           'KeyPressedCallback',    @(h,ev)SliderCallback(h, ev, 'ResectZ'));
                panelResect.add('hfill', jSliderResectZ);   
            jPanelSurfaceResect.add('hfill', panelResect);
            
            % === HEMISPHERES SELECTION ===
            jToggleResectLeft  = gui_component('toggle', jPanelSurfaceResect, 'br center', 'Left',  {Insets(0,0,0,0), Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT)}, '', @ButtonResectLeftToggle_Callback);           
            jToggleResectRight = gui_component('toggle', jPanelSurfaceResect, '',          'Right', {Insets(0,0,0,0), Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT)}, '', @ButtonResectRightToggle_Callback);           
            gui_component('button', jPanelSurfaceResect, '',          'Reset', {Insets(0,0,0,0), Dimension(BUTTON_WIDTH, DEFAULT_HEIGHT)}, '', @ButtonResectResetCallback);
        jPanelOptions.add('br hfill', jPanelSurfaceResect);
 
        % Labels
        gui_component('label', jPanelOptions, 'br', '    Vertices : ');
        jLabelNbVertices = gui_component('label', jPanelOptions, '', '0');
        gui_component('label', jPanelOptions, '', '    Faces : ');
        jLabelNbFaces = gui_component('label', jPanelOptions, '', '0');

    jScrollPaneOptions = JScrollPane(jPanelOptions);
    jScrollPaneOptions.setBorder([]);
    jPanelNew.add(jScrollPaneOptions, BorderLayout.CENTER);
    
    % Create the BstPanel object that is returned by the function
    % => constructor BstPanel(jHandle, panelName, sControls)
    bstPanelNew = BstPanel(panelName, ...
                           jPanelNew, ...
                           struct('jToolbar',       jToolbar, ...
                                  'jPanelOptions',          jPanelOptions, ...
                                  'jPanelSurfaceOptions',   jPanelSurfaceOptions, ...
                                  'jPanelSurfaceResect',    jPanelSurfaceResect, ...
                                  'jPanelDataOptions',      jPanelDataOptions, ...                             
                                  'jLabelNbVertices',       jLabelNbVertices, ...
                                  'jLabelNbFaces',          jLabelNbFaces, ...
                                  'jSliderSurfAlpha',       jSliderSurfAlpha, ...
                                  'jLabelSurfAlpha',        jLabelSurfAlpha, ...
                                  'jButtonSurfColor',       jButtonSurfColor, ...
                                  'jLabelSurfSmoothValue',  jLabelSurfSmoothValue, ...
                                  'jSliderSurfSmoothValue', jSliderSurfSmoothValue, ...
                                  'jButtonSurfSulci',       jButtonSurfSulci, ...
                                  'jButtonSurfEdge',        jButtonSurfEdge, ...
                                  'jSliderResectX',         jSliderResectX, ...
                                  'jSliderResectY',         jSliderResectY, ...
                                  'jSliderResectZ',         jSliderResectZ, ...
                                  'jToggleResectLeft',      jToggleResectLeft, ...
                                  'jToggleResectRight',     jToggleResectRight, ...
                                  'jLabelAlphaTitle',       jLabelAlphaTitle, ...
                                  'jLabelDataAlphaTitle',   jLabelDataAlphaTitle, ...
                                  'jSliderDataAlpha',       jSliderDataAlpha, ...
                                  'jLabelDataAlpha',        jLabelDataAlpha, ...
                                  'jLabelSizeTitle',        jLabelSizeTitle, ...
                                  'jLabelSize',             jLabelSize, ...
                                  'jSliderSize',            jSliderSize, ...
                                  'sliderSizeVector',       sliderSizeVector, ...
                                  'jLabelThreshTitle',      jLabelThreshTitle, ...
                                  'jSliderDataThresh',      jSliderDataThresh, ...
                                  'jLabelDataThresh',       jLabelDataThresh));


    %% ===== SLIDER QUICK PREVIEW =====
    function SliderQuickPreview(jSlider, jText, isPercent)
        if (jSlider == jSliderSize)
            jText.setText(sprintf('%d', sliderSizeVector(double(jSlider.getValue()))));
        elseif isPercent
            jText.setText(sprintf('%d%%', double(jSlider.getValue())));
        else
            jText.setText(sprintf('%d', double(jSlider.getValue())));
        end
    end
                              
    %% ===== RESET RESECT CALLBACK =====
    function ButtonResectResetCallback(varargin)
        import java.awt.event.*;
        % Reset initial resect sliders positions
        jSliderResectX.setValue(0);
        jSliderResectY.setValue(0);
        jSliderResectZ.setValue(0);
        SliderCallback([], MouseEvent(jSliderResectX, 0, 0, 0, 0, 0, 1, 0), 'ResectX');
        SliderCallback([], MouseEvent(jSliderResectY, 0, 0, 0, 0, 0, 1, 0), 'ResectY');
        SliderCallback([], MouseEvent(jSliderResectZ, 0, 0, 0, 0, 0, 1, 0), 'ResectZ');
    end

    %% ===== RESECT LEFT TOGGLE CALLBACK =====
    function ButtonResectLeftToggle_Callback(varargin)
        if jToggleResectLeft.isSelected()
            jToggleResectRight.setSelected(0);
            SelectHemispheres('left');
        else
            SelectHemispheres('none');
        end
    end

    %% ===== RESECT RIGHT TOGGLE CALLBACK =====
    function ButtonResectRightToggle_Callback(varargin)
        if jToggleResectRight.isSelected()
            jToggleResectLeft.setSelected(0);
            SelectHemispheres('right');
        else
            SelectHemispheres('none');
        end
    end
end


%% =================================================================================
%  === CONTROLS CALLBACKS  =========================================================
%  =================================================================================
%% ===== SLIDERS CALLBACKS =====
function SliderCallback(hObject, event, target)
    % Get panel controls
    ctrl = bst_get('PanelControls', 'Surface');
    % Get slider pointer
    jSlider = event.getSource();
    % If slider is not enabled : do nothing
    if ~jSlider.isEnabled()
        return
    end

    % Get handle to current 3DViz figure
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get current surface index (in the current figure)
    iSurface = getappdata(hFig, 'iSurface');
    % If surface data is not accessible
    if isempty(hFig) || isempty(iSurface)
        return;
    end
    % Get figure AppData (figure's surfaces configuration)
    TessInfo = getappdata(hFig, 'Surface');
    if (iSurface > length(TessInfo))
        return;
    end
    % Is selected surface a MRI/slices surface
    isAnatomy = strcmpi(TessInfo(iSurface).Name, 'Anatomy');
    
    % Get slider value and update surface value
    switch (target)
        case 'SurfAlpha'
            % Update value in Surface array
            TessInfo(iSurface).SurfAlpha = jSlider.getValue() / 100;
            % Display value in the label associated with the slider
            ctrl.jLabelSurfAlpha.setText(sprintf('%d%%', round(TessInfo(iSurface).SurfAlpha * 100)));
            % Update current surface
            setappdata(hFig, 'Surface', TessInfo);
            % For MRI: redraw all slices
            if isAnatomy
                figure_callback(hFig, 'UpdateMriDisplay', hFig, [], TessInfo, iSurface);
            % Else: Update color display on the surface
            else
                figure_callback(hFig, 'UpdateSurfaceAlpha', hFig, iSurface);
            end
    
        case 'SurfSmoothValue'
            SurfSmoothValue = jSlider.getValue() / 100;
            SetSurfaceSmooth(hFig, iSurface, SurfSmoothValue);

        case 'DataAlpha'
            % Update value in Surface array
            TessInfo(iSurface).DataAlpha = jSlider.getValue() / 100;
            ctrl.jLabelDataAlpha.setText(sprintf('%d%%', round(TessInfo(iSurface).DataAlpha * 100)));
            % Update current surface
            setappdata(hFig, 'Surface', TessInfo);
            % Update color display on the surface
            figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurface);
            % Set the new value as he default value (NOT FOR MRI)
            DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
            DefaultSurfaceDisplay.DataAlpha = TessInfo(iSurface).DataAlpha;
            bst_set('DefaultSurfaceDisplay', DefaultSurfaceDisplay); 
            
        case 'DataThreshold'
            % Update value in Surface array
            TessInfo(iSurface).DataThreshold = jSlider.getValue() / 100;
            ctrl.jLabelDataThresh.setText(sprintf('%d%%', round(TessInfo(iSurface).DataThreshold * 100)));
            % Update current surface
            setappdata(hFig, 'Surface', TessInfo);
            % Update color display on the surface
            figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurface);
            % Set the new value as the default value (NOT FOR MRI)
            DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
            DefaultSurfaceDisplay.DataThreshold = TessInfo(iSurface).DataThreshold;
            bst_set('DefaultSurfaceDisplay', DefaultSurfaceDisplay); 
            
        case 'SizeThreshold'
            % Update value in Surface array
            TessInfo(iSurface).SizeThreshold = ctrl.sliderSizeVector(jSlider.getValue());
            ctrl.jLabelSize.setText(sprintf('%d', round(TessInfo(iSurface).SizeThreshold)));
            % Update current surface
            setappdata(hFig, 'Surface', TessInfo);
            % Update color display on the surface
            figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurface);
            % Set the new value as the default value (NOT FOR MRI)
            if ~isAnatomy
                DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
                DefaultSurfaceDisplay.SizeThreshold = TessInfo(iSurface).SizeThreshold;
                bst_set('DefaultSurfaceDisplay', DefaultSurfaceDisplay); 
            end
        case {'ResectX', 'ResectY', 'ResectZ'}
            % Get target axis
            dim = find(strcmpi(target, {'ResectX', 'ResectY', 'ResectZ'}));
            % JSliderResect values : [-100,100]
            if isAnatomy
                % Get MRI size
                sMri = bst_memory('GetMri', TessInfo(iSurface).SurfaceFile);
                cubeSize = size(sMri.Cube);
                % Change slice position
                newPos = round((jSlider.getValue()+100) / 200 * cubeSize(dim));
                newPos = bst_saturate(newPos, [1, cubeSize(dim)]);
                TessInfo(iSurface).CutsPosition(dim) = newPos;
                % Update MRI display
                figure_callback(hFig, 'UpdateMriDisplay', hFig, dim, TessInfo, iSurface);
            else
                ResectSurface(hFig, iSurface, dim, jSlider.getValue() / 100);
            end

        otherwise
            error('Unknow slider');
    end
end



%% ===== BUTTON SURFACE COLOR CALLBACK =====
function ButtonSurfColorCallback(varargin)
    % Get handle to current 3DViz figure
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get figure AppData (figure's surfaces configuration)
    TessInfo = getappdata(hFig, 'Surface');
    % Get current surface index (in the current figure)
    iSurface = getappdata(hFig, 'iSurface');
    % Ignore MRI slices
    if strcmpi(TessInfo(iSurface).Name, 'Anatomy')
        return
    end
    % Ask user to select a color
    colorCortex = uisetcolor(TessInfo(iSurface).AnatomyColor(2,:), 'Select surface color');
    if (length(colorCortex) ~= 3)
        return
    end
    % Change surface color
    SetSurfaceColor(hFig, iSurface, colorCortex);
end
             


%% ===== BUTTON "SULCI" CALLBACK =====
function ButtonShowSulciCallback(hObject, event)
    % Get handle to current 3DViz figure
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get current surface index (in the current figure)
    iSurface = getappdata(hFig, 'iSurface');
    % Get handle to "View" button
    jButtonSurfSulci = event.getSource();
    % Show/hide sulci map in figure display
    SetShowSulci(hFig, iSurface, jButtonSurfSulci.isSelected());
    % Set the new value as the default value
    DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
    DefaultSurfaceDisplay.SurfShowSulci = jButtonSurfSulci.isSelected();
    bst_set('DefaultSurfaceDisplay', DefaultSurfaceDisplay);
end

%% ===== SET SHOW SULCI =====
% Usage : SetShowSulci(hFig, iSurfaces, status)
% Parameters : 
%     - hFig : handle to a 3DViz figure
%     - iSurfaces : can be a single indice or an array of indices
%     - status    : 1=display, 0=hide
function SetShowSulci(hFig, iSurfaces, status)
    % Get surfaces list 
    TessInfo = getappdata(hFig, 'Surface');
    % Process all surfaces
    for iSurf = iSurfaces
        % Shet status : show/hide
        TessInfo(iSurf).SurfShowSulci = status;
    end
    % Update figure's AppData (surfaces configuration)
    setappdata(hFig, 'Surface', TessInfo);
    % Update surface display
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurf);
end


%% ===== SHOW SURFACE EDGES =====
function ButtonShowEdgesCallback(varargin)
    % Get handle to current 3DViz figure
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get current surface (in the current figure)
    TessInfo = getappdata(hFig, 'Surface');
    iSurf    = getappdata(hFig, 'iSurface');
    % Set edges display on/off
    TessInfo(iSurf).SurfShowEdges = ~TessInfo(iSurf).SurfShowEdges;
    setappdata(hFig, 'Surface', TessInfo);
    % Update display
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurf);
end


%% ===== HEMISPHERE SELECTION RADIO CALLBACKS =====
function SelectHemispheres(name)
    % Get handle to current 3DViz figure
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get surface properties
    TessInfo = getappdata(hFig, 'Surface');
    iSurf    = getappdata(hFig, 'iSurface');
    % Ignore MRI
    if strcmpi(TessInfo(iSurf).Name, 'Anatomy')
        return;
    end
    % Update surface Resect field
    TessInfo(iSurf).Resect = name;
    setappdata(hFig, 'Surface', TessInfo);
    
    % Reset all the resect sliders
    ctrl = bst_get('PanelControls', 'Surface');
    ctrl.jSliderResectX.setValue(0);
    ctrl.jSliderResectY.setValue(0);
    ctrl.jSliderResectZ.setValue(0);
    
    % Display progress bar
    bst_progress('start', 'Select hemisphere', 'Selecting hemisphere...');
    % Update surface display
    figure_callback(hFig, 'UpdateSurfaceAlpha', hFig, iSurf);
    % Display progress bar
    bst_progress('stop');
end


%% ===== RESECT SURFACE =====
function ResectSurface(hFig, iSurf, resectDim, resectValue)
    % Get surfaces description
    TessInfo = getappdata(hFig, 'Surface');
    % If previously using "Select hemispheres"
    if ischar(TessInfo(iSurf).Resect)
        % Reset "Resect" field
        TessInfo(iSurf).Resect = [0 0 0];
    end
    % Update value in Surface array
    TessInfo(iSurf).Resect(resectDim) = resectValue;
    % Update surface
    setappdata(hFig, 'Surface', TessInfo);
    % Hide trimmed part of the surface
    figure_callback(hFig, 'UpdateSurfaceAlpha', hFig, iSurf);
    
    % Deselect both Left and Right buttons
    ctrl = bst_get('PanelControls', 'Surface');
    ctrl.jToggleResectLeft.setSelected(0);
    ctrl.jToggleResectRight.setSelected(0);
end


%% ===== ADD SURFACE CALLBACK =====
function ButtonAddSurfaceCallback(varargin)
    % Get target figure handle
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get Current subject
    SubjectFile = getappdata(hFig, 'SubjectFile');
    if isempty(SubjectFile)
        return
    end
    sSubject = bst_get('Subject', SubjectFile);
    if isempty(sSubject)
        return
    end

    % List of available surfaces types
    typesList = {};
    if ~isempty(sSubject.iScalp)
        typesList{end+1} = 'Scalp';
    end
    if ~isempty(sSubject.iOuterSkull)
        typesList{end+1} = 'OuterSkull';
    end
    if ~isempty(sSubject.iInnerSkull)
        typesList{end+1} = 'InnerSkull';
    end
    if ~isempty(sSubject.iCortex)
        typesList{end+1} = 'Cortex';
    end
    if ~isempty(sSubject.iAnatomy)
        typesList{end+1} = 'Anatomy';
    end
    if isempty(typesList)
        return
    end
    % Ask user which kind of surface he wants to add to the figure 3DViz
    surfaceType = java_dialog('question', 'What kind of surface would you like to display ?', 'Add surface', [], typesList, typesList{1});

    % Switch between surfaces types
    switch (surfaceType)
        case 'Anatomy'
            SurfaceFile = sSubject.Anatomy(sSubject.iAnatomy(1)).FileName;
        case 'Cortex'
            SurfaceFile = sSubject.Surface(sSubject.iCortex(1)).FileName;
        case 'Scalp'
            SurfaceFile = sSubject.Surface(sSubject.iScalp(1)).FileName;
        case 'InnerSkull'
            SurfaceFile = sSubject.Surface(sSubject.iInnerSkull(1)).FileName;
        case 'OuterSkull'
            SurfaceFile = sSubject.Surface(sSubject.iOuterSkull(1)).FileName;
        otherwise
            return;
    end
    % Get surface list
    TessInfo = getappdata(hFig, 'Surface');
    % Add surface to the figure
    iTess = AddSurface(hFig, SurfaceFile);
    % 3D MRI: Update Colormap
    if strcmpi(surfaceType, 'Anatomy')
        % Get figure
        [hFig,iFig,iDS] = bst_figures('GetFigure', hFig);
        % Update colormap
        figure_3d('ColormapChangedCallback', iDS, iFig);
    end
    % Reload scouts (only if new surface was added)
    if (iTess > length(TessInfo))
        panel_scouts('ReloadScouts', hFig);
    end
end


%% ===== REMOVE SURFACE CALLBACK =====
function ButtonRemoveSurfaceCallback(varargin)
    % Get target figure handle
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get current surface index
    iSurface = getappdata(hFig, 'iSurface');
    if isempty(iSurface)
        return
    end
    % Remove surface
    RemoveSurface(hFig, iSurface);
    % Update "Surfaces" panel
    UpdatePanel();
end


%% ===== SURFACE BUTTON CLICKED CALLBACK =====
function ButtonSurfaceClickedCallback(hObject, event, varargin)
    % Get current 3DViz figure
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get index of the surface associated to this button
    iSurface = str2num(event.getSource.getName());
    % Store current surface index 
    setappdata(hFig, 'iSurface', iSurface);
    % Update surface properties
    UpdateSurfaceProperties();
    % Reload scouts
    panel_scouts('ReloadScouts', hFig);
end



%% =================================================================================
%  === EXTERNAL CALLBACKS  =========================================================
%  =================================================================================
%% ===== UPDATE PANEL =====
function UpdatePanel(varargin)
    % Get panel controls
    panelSurfacesCtrl = bst_get('PanelControls', 'Surface');
    if isempty(panelSurfacesCtrl)
        return
    end
    % If no current 3D figure defined
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        % Remove surface buttons
        CreateSurfaceList(panelSurfacesCtrl.jToolbar, 0);
        % Disable all panel controls
        gui_enable([panelSurfacesCtrl.jToolbar, panelSurfacesCtrl.jPanelOptions], 0);  
    else
        % Enable Surfaces selection panel
        gui_enable(panelSurfacesCtrl.jToolbar, 1);
        % Update surfaces list
        nbSurfaces = CreateSurfaceList(panelSurfacesCtrl.jToolbar, hFig);
        % If no surface is available
        if (nbSurfaces <= 0)
            % Disable "Display" and "Options" panel
            gui_enable(panelSurfacesCtrl.jPanelOptions, 0);
            % Else : one or more surfaces are available
        else
            % Enable "Display" and "Options" panel
            gui_enable(panelSurfacesCtrl.jPanelOptions, 1);
            % Update surface properties
            UpdateSurfaceProperties();
        end
    end
end

%% ===== CURRENT FREQ CHANGED CALLBACK =====
function CurrentFreqChangedCallback(iDS, iFig) %#ok<DEFNU>
    global GlobalData;
    % Get figure appdata
    hFig = GlobalData.DataSet(iDS).Figure(iFig).hFigure;
    TfInfo = getappdata(hFig, 'Timefreq');
    % If no frequencies in this figure
    if getappdata(hFig, 'isStaticFreq')
        return;
    end
    % If there are some time-frequency recordings in this
    if ~isempty(TfInfo)
        % Update frequency to display
        TfInfo.iFreqs = GlobalData.UserFrequencies.iCurrentFreq;
        setappdata(hFig, 'Timefreq', TfInfo);
        % Update display
        UpdateSurfaceData(hFig);
    end
end


%% ===== DISPATCH FIGURE CALLBACKS =====
function figure_callback(hFig, CallbackName, varargin)
    % Get figure type
    FigureId = getappdata(hFig, 'FigureId');
    % Different figure types
    switch (FigureId.Type)
        case 'MriViewer'
            figure_mri(CallbackName, varargin{:});
        case '3DViz'
            figure_3d(CallbackName, varargin{:});
    end
end


%% ===== CURRENT FIGURE CHANGED =====
function CurrentFigureChanged_Callback() %#ok<DEFNU>
    UpdatePanel();
end


%% ===== CREATE SURFACES LIST =====
function nbSurfaces = CreateSurfaceList(jToolbar, hFig)
    % Java initializations
    import java.awt.*;
    import javax.swing.*;
    import org.brainstorm.icon.*;

    nbSurfaces = 0;
    % Remove all toolbar surface buttons
    for iComp = 1:jToolbar.getComponentCount()-5
        jToolbar.remove(1);
    end
    % If no figure is specified : return
    if isempty(hFig) || ~ishandle(hFig) || (hFig == 0)
        return;
    end
    % Create a button group for Surfaces and "Add" button
    jButtonGroup = ButtonGroup();
    
    % If a figure is defined 
    if ishandle(hFig)
        % Get selected surface index
        iSurface = getappdata(hFig, 'iSurface');
        % Loop on all the available surfaces for this figure
        TessInfo = getappdata(hFig, 'Surface');
        for iSurf = 1:length(TessInfo)
            % Select only one button
            isSelected = (iSurf == iSurface);
            % Get button icon (depends on surface name)
            switch lower(TessInfo(iSurf).Name)
                case 'cortex'
                    iconButton = IconLoader.ICON_SURFACE_CORTEX;
                case 'scalp'
                    iconButton = IconLoader.ICON_SURFACE_SCALP;
                case 'innerskull'
                    iconButton = IconLoader.ICON_SURFACE_INNERSKULL;
                case 'outerskull'
                    iconButton = IconLoader.ICON_SURFACE_OUTERSKULL;
                case 'other'
                    iconButton = IconLoader.ICON_SURFACE;
                case 'anatomy'
                    iconButton = IconLoader.ICON_ANATOMY;
            end
            % Create surface button 
            jButtonSurf = JToggleButton(iconButton, isSelected);
            jButtonSurf.setMaximumSize(Dimension(24,24));
            jButtonSurf.setPreferredSize(Dimension(24,24));
            jButtonSurf.setToolTipText(TessInfo(iSurf).SurfaceFile);
            % Store the surface index as the button Name
            jButtonSurf.setName(sprintf('%d', iSurf));
            % Attach a click callback
            java_setcb(jButtonSurf, 'ActionPerformedCallback', @ButtonSurfaceClickedCallback);
            % Add button to button group
            jButtonGroup.add(jButtonSurf);
            % Add button to toolbar, at the end of the surfaces list
            iButton = jToolbar.getComponentCount() - 4;
            jToolbar.add(jButtonSurf, iButton);
        end
        % Return number of surfaces added
        nbSurfaces = length(TessInfo);
    else
        % No surface available for current figure
        nbSurfaces = 0;
    end
   
    % Update graphical composition of panel
    jToolbar.updateUI();
end


%% ===== UPDATE SURFACE PROPERTIES =====
function UpdateSurfaceProperties()
% disp('=== panel_surface > UpdateSurfaceProperties ===');
    import org.brainstorm.list.*;
    % Get current figure handle
    hFig = bst_figures('GetCurrentFigure', '3D');
    if isempty(hFig)
        return
    end
    % Get panel controls
    ctrl = bst_get('PanelControls', 'Surface');
    if isempty(ctrl)
        return
    end
    % Get selected surface properties
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo)
        return;
    end
    % Get selected surface index
    iSurface = getappdata(hFig, 'iSurface');
    % If surface is sliced MRI
    isAnatomy = strcmpi(TessInfo(iSurface).Name, 'Anatomy');

    % ==== Surface properties ====
    % Number of vertices
    ctrl.jLabelNbVertices.setText(sprintf('%d', TessInfo(iSurface).nVertices));
    % Number of faces
    ctrl.jLabelNbFaces.setText(sprintf('%d', TessInfo(iSurface).nFaces));
    % Surface alpha
    ctrl.jSliderSurfAlpha.setValue(100 * TessInfo(iSurface).SurfAlpha);
    ctrl.jLabelSurfAlpha.setText(sprintf('%d%%', round(100 * TessInfo(iSurface).SurfAlpha)));
    % Surface color
    surfColor = TessInfo(iSurface).AnatomyColor(2, :);
    ctrl.jButtonSurfColor.setBackground(java.awt.Color(surfColor(1),surfColor(2),surfColor(3)));
    % Surface smoothing ALPHA
    ctrl.jSliderSurfSmoothValue.setValue(100 * TessInfo(iSurface).SurfSmoothValue);
    ctrl.jLabelSurfSmoothValue.setText(sprintf('%d%%', round(100 * TessInfo(iSurface).SurfSmoothValue)));
    % Show sulci button
    ctrl.jButtonSurfSulci.setSelected(TessInfo(iSurface).SurfShowSulci);
    % Show surface edges button
    ctrl.jButtonSurfEdge.setSelected(TessInfo(iSurface).SurfShowEdges);
    
    % ==== Resect properties ====
    % Ignore for MRI slices
    if isAnatomy
        sMri = bst_memory('GetMri', TessInfo(iSurface).SurfaceFile);
        ResectXYZ = double(TessInfo(iSurface).CutsPosition) ./ size(sMri.Cube) * 200 - 100;
        radioSelected = 'none';
    elseif ischar(TessInfo(iSurface).Resect)
        ResectXYZ = [0,0,0];
        radioSelected = TessInfo(iSurface).Resect;
    else
        ResectXYZ = 100 * TessInfo(iSurface).Resect;
        radioSelected = 'none';
    end
    % X, Y, Z
    ctrl.jSliderResectX.setValue(ResectXYZ(1));
    ctrl.jSliderResectY.setValue(ResectXYZ(2));
    ctrl.jSliderResectZ.setValue(ResectXYZ(3));
    
    % Select one radio button
    switch (radioSelected)
        case 'left'
            ctrl.jToggleResectLeft.setSelected(1);
            ctrl.jToggleResectRight.setSelected(0);
        case 'right'
            ctrl.jToggleResectRight.setSelected(1);
            ctrl.jToggleResectLeft.setSelected(0);
        case 'none'
            ctrl.jToggleResectLeft.setSelected(0);
            ctrl.jToggleResectRight.setSelected(0);
    end
    
    % ==== Data properties ====
    % Enable/disable controls
    isOverlay = ~isempty(TessInfo(iSurface).DataSource.FileName);
    isOverlayStat = isOverlay && ismember(file_gettype(TessInfo(iSurface).DataSource.FileName), {'presults', 'pdata', 'ptimefreq'});
    gui_enable([ctrl.jLabelDataAlphaTitle, ctrl.jSliderDataAlpha, ctrl.jLabelDataAlpha, ...
                ctrl.jLabelSizeTitle, ctrl.jLabelSize, ctrl.jSliderSize], isOverlay, 0);
    gui_enable([ctrl.jLabelThreshTitle, ctrl.jSliderDataThresh, ctrl.jLabelDataThresh], isOverlay && ~isOverlayStat, 0);
    % Data threshold
    ctrl.jSliderDataThresh.setValue(100 * TessInfo(iSurface).DataThreshold);
    ctrl.jLabelDataThresh.setText(sprintf('%d%%', round(100 * TessInfo(iSurface).DataThreshold)));
    % Size threshold
    iSlider = be_closest(ctrl.sliderSizeVector, TessInfo(iSurface).SizeThreshold);
    ctrl.jSliderSize.setValue(iSlider);
    ctrl.jLabelSize.setText(sprintf('%d', TessInfo(iSurface).SizeThreshold));
    % Data alpha
    ctrl.jSliderDataAlpha.setValue(100 * TessInfo(iSurface).DataAlpha);
    ctrl.jLabelDataAlpha.setText(sprintf('%d%%', round(100 * TessInfo(iSurface).DataAlpha)));
end


%% ===== ADD A SURFACE =====
% Add a surface to a given 3DViz figure
% USAGE : iTess = panel_surface('AddSurface', hFig, surfaceFile)
% OUTPUT: Indice of the surface in the figure's surface array
function iTess = AddSurface(hFig, surfaceFile)
    % ===== CHECK EXISTENCE =====
    % Check whether filename is an absolute or relative path
    surfaceFile = file_short(surfaceFile);
    % Get figure appdata (surfaces configuration)
    TessInfo = getappdata(hFig, 'Surface');
    % Check that this surface is not already displayed in 3DViz figure
    iTess = find(file_compare({TessInfo.SurfaceFile}, surfaceFile));
    if ~isempty(iTess)
        disp('BST> This surface is already displayed. Ignoring...');
        return
    end
    % Get figure type
    FigureId = getappdata(hFig, 'FigureId');
    % Progress bar
    isNewProgressBar = ~bst_progress('isVisible');
    bst_progress('start', 'Add surface', 'Updating display...');
    
    % ===== BUILD STRUCTURE =====
    % Add a new surface at the end of the figure's surfaces list
    iTess = length(TessInfo) + 1;
    TessInfo(iTess) = db_template('TessInfo');                       
    % Set the surface properties
    TessInfo(iTess).SurfaceFile = surfaceFile;
    TessInfo(iTess).DataSource.Type     = '';
    TessInfo(iTess).DataSource.FileName = '';

    % ===== PLOT OBJECT =====
    % Get file type (tessalation or MRI)
    fileType = file_gettype(surfaceFile);
    % === TESSELATION ===
    if any(strcmpi(fileType, {'cortex','scalp','innerskull','outerskull','tess'}))
        % === LOAD SURFACE ===
        % Load surface file
        sSurface = bst_memory('LoadSurface', surfaceFile);
        % Get some properties
        TessInfo(iTess).Name      = sSurface.Name;
        TessInfo(iTess).nVertices = size(sSurface.Vertices, 1);
        TessInfo(iTess).nFaces    = size(sSurface.Faces, 1);

        % === PLOT SURFACE ===
        switch (FigureId.Type)
            case 'MriViewer'
                % Nothing to do: surface will be displayed as an overlay slice in figure_mri.m
            case {'3DViz', 'Topography'}
                % Create and display surface patch
                [hFig, TessInfo(iTess).hPatch] = figure_3d('PlotSurface', hFig, ...
                                         sSurface.Faces, ...
                                         sSurface.Vertices, ...
                                         TessInfo(iTess).AnatomyColor(2,:), ...
                                         TessInfo(iTess).SurfAlpha);
        end
        % Update figure's surfaces list and current surface pointer
        setappdata(hFig, 'Surface',  TessInfo);
        
        % Show sulci map if needed 
        if TessInfo(iTess).SurfShowSulci
            SetShowSulci(hFig, iTess, 1);
        end
        % If displaying the first surface, and it is a cortex: unzoom a bit
        if (iTess == 1)
            if strcmpi(sSurface.Name, 'Cortex')
                zoom(hFig, 0.87);
                zoom reset;
            elseif strcmpi(sSurface.Name, 'Scalp')
                zoom(hFig, 1.1);
                zoom reset;
            end
        end
        
    % === MRI ===
    elseif strcmpi(fileType, 'subjectimage')
        % === LOAD MRI ===
        sMri = bst_memory('LoadMri', surfaceFile);
        if isempty(sMri)
            return
        end
        TessInfo(iTess).Name = 'Anatomy';
        % Initial position of the cuts : middle in each direction
        TessInfo(iTess).CutsPosition = round(size(sMri.Cube) / 2);
        TessInfo(iTess).SurfSmoothValue = .3;
        % Colormap
        TessInfo(iTess).ColormapType = 'anatomy';
        bst_colormaps('AddColormapToFigure', hFig, TessInfo(iTess).ColormapType);
        % Update figure's surfaces list and current surface pointer
        setappdata(hFig, 'Surface',  TessInfo);

        % === PLOT MRI ===
        switch (FigureId.Type)
            case 'MriViewer'
                % Configure MRIViewer
                figure_mri('SetupMri', hFig);
            case '3DViz'
                % Camera basic orientation: TOP
                figure_3d('SetStandardView', hFig, 'top');
        end
        % Plot MRI
        PlotMri(hFig);
    end
    % Update default surface
    setappdata(hFig, 'iSurface', iTess);
    % Automatically set transparencies (to view different layers at the same time)
    SetAutoTransparency(hFig);
    % Close progress bar
    drawnow;
    if isNewProgressBar
        bst_progress('stop');
    end
    % Update panel
    UpdatePanel();
end
   


%% ===== SET DATA SOURCE FOR A SURFACE =====
%Associate a data/results matrix to a surface.
% Usage : SetSurfaceData(hFig, iTess, dataType, dataFile, isStat)
% Parameters : 
%     - hFig : handle to a 3DViz figure
%     - iTess        : indice of the surface to update (in hFig appdata)
%     - dataType     : type of data to overlay on the surface {'Source', 'Data', ...}
%     - dataFile     : filename of the data to display over the surface
%     - isStat       : 1, if results is a statistical result; 0, else
function isOk = SetSurfaceData(hFig, iTess, dataType, dataFile, isStat) %#ok<DEFNU>
    global GlobalData;
    % Get figure index in DataSet figures list
    [tmp__, iFig, iDS] = bst_figures('GetFigure', hFig);
    if isempty(iDS)
        error('No DataSet acessible for this 3D figure');
    end
    % Get surfaces list for this figure
    TessInfo = getappdata(hFig, 'Surface');
    isAnatomy = strcmpi(TessInfo(iTess).Name, 'Anatomy');
    
    % === GET DATA THRESHOLD ===
    % Get defaults for surface display
    DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
    % Cortex
    %if strcmpi(TessInfo(iTess).Name, 'Cortex') && ~isStat
    if ~isStat
        % Data/size threshold
        dataThreshold = DefaultSurfaceDisplay.DataThreshold;
        if isAnatomy
            sizeThreshold = 1;
        else
            sizeThreshold = DefaultSurfaceDisplay.SizeThreshold;
        end
    % Anatomy or Statistics : 0%
    elseif isAnatomy || isStat
        dataThreshold = 0;
        sizeThreshold = 1;
    % Else: normal data on scalp
    else
        dataThreshold = 0.5;
        sizeThreshold = 1;
    end
    % Static figure
    setappdata(hFig, 'isStatic', (GlobalData.DataSet(iDS).Measures.NumberOfSamples <= 2));
    
    % === PREPARE SURFACE ===
    TessInfo(iTess).DataSource.Type     = dataType;
    TessInfo(iTess).DataSource.FileName = dataFile;
    TessInfo(iTess).DataThreshold       = dataThreshold;
    TessInfo(iTess).SizeThreshold       = sizeThreshold;
    TessInfo(iTess).DataAlpha           = DefaultSurfaceDisplay.DataAlpha;
    % Type of data displayed on the surface: sources/recordings/nothing
    switch (dataType)
        case 'Data'
            % Get loaded data
            iDS = bst_memory('GetDataSetData', dataFile);
            % Select appropriate colormap
            if ~isempty(GlobalData.DataSet(iDS).Measures.ColormapType)
                ColormapType = GlobalData.DataSet(iDS).Measures.Colormap;
            elseif isStat
                ColormapType = 'stat2';
            else
                ColormapType = 'eeg';
            end
            setappdata(hFig, 'DataFile', dataFile);
            
        case 'Source'
            % Get loaded results
            [iDS, iRes] = bst_memory('GetDataSetResult', dataFile);
            % Select appropriate colormap
            if ~isempty(GlobalData.DataSet(iDS).Results(iRes).ColormapType)
                ColormapType = GlobalData.DataSet(iDS).Results(iRes).ColormapType;
            elseif isStat
                ColormapType = 'stat1';
            else
                ColormapType = 'source';
            end
            % Copy the surface atlas
            TessInfo(iTess).DataSource.Atlas = GlobalData.DataSet(iDS).Results(iRes).Atlas;
            setappdata(hFig, 'ResultsFile', dataFile);
            
        case 'Dipoles'
            ColormapType = 'source';
            panel_dipoles('AddDipoles', hFig, dataFile, 0);
            
        case 'Timefreq'
            % Get study
            [sStudy, iStudy, iTf] = bst_get('TimefreqFile', dataFile);
            if isempty(sStudy)
                error('File is not registered in database.');
            end
            % Get loaded time-freq structure
            [iDS, iTimefreq] = bst_memory('LoadTimefreqFile', sStudy.Timefreq(iTf).FileName);
             % Set "Static" status for this figure
            setappdata(hFig, 'isStatic', (GlobalData.DataSet(iDS).Timefreq(iTimefreq).NumberOfSamples <= 2));
            isStaticFreq = (size(GlobalData.DataSet(iDS).Timefreq(iTimefreq).TF,3) <= 1);
            setappdata(hFig, 'isStaticFreq', isStaticFreq);

            % Create options structure
            TfInfo = db_template('TfInfo');
            TfInfo.FileName = sStudy.Timefreq(iTf).FileName;
            TfInfo.Comment  = sStudy.Timefreq(iTf).Comment;
            TfInfo.RowName  = [];
            if isStaticFreq
                TfInfo.iFreqs = [];
            elseif ~isempty(GlobalData.UserFrequencies.iCurrentFreq)
                TfInfo.iFreqs = GlobalData.UserFrequencies.iCurrentFreq;
            else
                TfInfo.iFreqs = 1;
            end
            % Data type
            switch (GlobalData.DataSet(iDS).Timefreq(iTimefreq).Method)
                case 'morlet',   TfInfo.Function = 'power';       ColormapType = 'timefreq';
                case 'fft',      TfInfo.Function = 'power';       ColormapType = 'timefreq';
                case 'psd',      TfInfo.Function = 'power';       ColormapType = 'timefreq';
                case 'hilbert',  TfInfo.Function = 'magnitude';   ColormapType = 'timefreq';
                case 'instfreq', TfInfo.Function = 'other';       ColormapType = 'timefreq';
                case 'corr',     TfInfo.Function = 'other';       ColormapType = 'connect1';
                case 'cohere',   TfInfo.Function = 'other';       ColormapType = 'connect1';
                case 'granger',  TfInfo.Function = 'other';       ColormapType = 'connect1';
                case 'plv',      TfInfo.Function = 'magnitude';   ColormapType = 'connect1';
                case 'plvt',     TfInfo.Function = 'magnitude';   ColormapType = 'connect1';
                otherwise,       TfInfo.Function = 'power';       ColormapType = 'timefreq';
            end
            % Set figure data
            setappdata(hFig, 'Timefreq', TfInfo);
            % Display options panel
            isDisplayTab = ~strcmpi(TfInfo.Function, 'other');
            if isDisplayTab
                gui_brainstorm('ShowToolTab', 'Display');
            end
            % Get atlas if it is a source file
            if strcmpi(GlobalData.DataSet(iDS).Timefreq(iTimefreq).DataType, 'results') && ~isempty(GlobalData.DataSet(iDS).Timefreq(iTimefreq).DataFile)
                % Get results index
                iResult = bst_memory('GetResultInDataSet', iDS, GlobalData.DataSet(iDS).Timefreq(iTimefreq).DataFile);
                % Get atlas
                if ~isempty(iResult) && ~isempty(GlobalData.DataSet(iDS).Results(iResult).Atlas) && ~isempty(GlobalData.DataSet(iDS).Results(iResult).Atlas.Scouts)
                    TessInfo(iTess).DataSource.Atlas = GlobalData.DataSet(iDS).Results(iResult).Atlas;
                end
                %setappdata(hFig, 'ResultsFile', GlobalData.DataSet(iDS).Timefreq(iTimefreq).DataFile);
            end

        case 'Surface'
            ColormapType = 'overlay';
            
        otherwise
            ColormapType = '';
            TessInfo(iTess).Data = [];
            TessInfo(iTess).DataWmat = [];
    end
    % Add colormap of the surface to the figure
    if ~isempty(ColormapType)
        TessInfo(iTess).ColormapType = ColormapType;
        bst_colormaps('AddColormapToFigure', hFig, ColormapType);
    end
    % Update figure appdata
    setappdata(hFig, 'Surface', TessInfo);
    % Plot surface
    isOk = UpdateSurfaceData(hFig, iTess);
    % Update  panel
    UpdatePanel();
end



%% ===== UPDATE SURFACE DATA =====
% Update the 'Data' field for given surfaces :
%    - Load data/results matrix (F, ImageGridAmp, ...) if it is not loaded yet
%    - Store global minimum/maximum of data
%    - Interpolate data matrix over the target surface if number of vertices does not match
%    - And update color display (ColormapChangedCallback)
%
% Usage:  UpdateSurfaceData(hFig, iSurfaces)
%         UpdateSurfaceData(hFig)
function isOk = UpdateSurfaceData(hFig, iSurfaces)
% disp('=== panel_surface > UpdateSurfaceData ===');
    global GlobalData;
    isOk = 1;
    % Get surfaces list 
    TessInfo = getappdata(hFig, 'Surface');
    % If the aim is to update all the surfaces 
    if (nargin < 2) || isempty(iSurfaces)
        iSurfaces = find(~cellfun(@(c)isempty(c.Type), {TessInfo.DataSource}));
        if isempty(iSurfaces)
            return
        end
    end
        
    % Get figure index (in DataSet structure)
    [tmp__, iFig, iDS] = bst_figures('GetFigure', hFig);
    % Find the DataSet indice that corresponds to the current figure
    if isempty(iDS)
        error('No DataSet acessible for this 3D figure');
    end
    
    % For each surface
    for iTess = iSurfaces
        % If surface patch object doesn't exist => error
        if isempty(TessInfo(iTess).hPatch)
            error('Patch is not displayed');
        end
        
        % ===== GET SURFACE DATA =====
        % Switch between different data types to display on the surface
        switch (TessInfo(iTess).DataSource.Type)
            case 'Data'
                % Get TimeVector and current time indice
                [TimeVector, CurrentTimeIndex] = bst_memory('GetTimeVector', iDS);
                % If surface is displayed : update it
                if ~isempty(TessInfo(iTess).hPatch) && ishandle(TessInfo(iTess).hPatch)
                    % Get vertices of surface
                    Vertices = get(TessInfo(iTess).hPatch, 'Vertices');
                    % Get selected channels indices and location
                    selChan = GlobalData.DataSet(iDS).Figure(iFig).SelectedChannels;
                    % Get sensors positions
                    chan_loc = figure_3d('GetChannelPositions', iDS, selChan);

                    % Interpolate data on scalp surface (only if Matrix is not computed yet, or channels changed)
                    % => TRICK : it is difficult to test if the sensors locations or the surface vertices changed
                    %            => Just the the number of channels and vertices (should be ok...)
                    if isempty(TessInfo(iTess).DataWmat) || ...
                            (size(TessInfo(iTess).DataWmat,2) ~= length(selChan)) || ...
                            (size(TessInfo(iTess).DataWmat,1) ~= length(Vertices))
                        switch lower(GlobalData.DataSet(iDS).Figure(iFig).Id.Modality)
                            case 'eeg',  excludeParam = .3;
                            case 'ecog', excludeParam = .3;
                            case 'seeg', excludeParam = .3;
                            case 'meg',  excludeParam = .5;
                            otherwise,   excludeParam = 0;
                        end
                        nbNeigh = 4;
                        TessInfo(iTess).DataWmat = bst_shepards(Vertices, chan_loc, nbNeigh, excludeParam);
                    end
                    % Set data for current time frame
                    % FT 11-Jan-10: Remove "single"
                    TessInfo(iTess).Data = TessInfo(iTess).DataWmat * GlobalData.DataSet(iDS).Measures.F(selChan, CurrentTimeIndex);
                    % Store minimum and maximum of displayed data
                    TessInfo(iTess).DataMinMax = [min(min(TessInfo(iTess).Data)), ...
                                                  max(max(TessInfo(iTess).Data))];
                end
                % Update "Static" status for this figure
                setappdata(hFig, 'isStatic', (GlobalData.DataSet(iDS).Measures.NumberOfSamples <= 2));

            case 'Source'
                % === LOAD RESULTS VALUES ===
                % Get results index
                iResult = bst_memory('GetResultInDataSet', iDS, TessInfo(iTess).DataSource.FileName);
                % If Results file is not found in GlobalData structure
                if isempty(iResult)
                    isOk = 0;
                    return
                end
                % If data matrix is not loaded for this file
                if isempty(GlobalData.DataSet(iDS).Results(iResult).ImageGridAmp) && isempty(GlobalData.DataSet(iDS).Results(iResult).ImagingKernel)
                    bst_memory('LoadResultsMatrix', iDS, iResult);
                end
                
                % === GET CURRENT VALUES ===
                % Get results values
                TessInfo(iTess).Data = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');
                % If min/max values for this file were not computed yet
                if isempty(TessInfo(iTess).DataMinMax)
                    TessInfo(iTess).DataMinMax = bst_memory('GetResultsMaximum', iDS, iResult);
                end
                % Reset Overlay cube
                TessInfo(iTess).OverlayCube = [];

                % Check the consistency between the number of results points (number of sources)
                % and the number of vertices of the target surface patch (IGNORE TEST FOR MRI)
                if strcmpi(TessInfo(iTess).Name, 'Anatomy')
                    % Nothing to check right now
                elseif ~isempty(TessInfo(iTess).DataSource.Atlas) 
                    if (length(TessInfo(iTess).Data) ~= length(TessInfo(iTess).DataSource.Atlas.Scouts))
                        bst_error(sprintf(['Number of sources (%d) is different from number of scouts (%d).\n\n' ...
                                  'Please compute the sources again.'], size(TessInfo(iTess).Data, 1), TessInfo(iTess).DataSource.Atlas.Scouts), 'Data mismatch', 0);
                        isOk = 0;
                        return;
                    end
                elseif (length(TessInfo(iTess).Data) ~= TessInfo(iTess).nVertices)
                    bst_error(sprintf(['Number of sources (%d) is different from number of vertices (%d).\n\n' ...
                              'Please compute the sources again.'], size(TessInfo(iTess).Data, 1), TessInfo(iTess).nVertices), 'Data mismatch', 0);
                    isOk = 0;
                    return;
                end
                % Update "Static" status for this figure
                setappdata(hFig, 'isStatic', (GlobalData.DataSet(iDS).Results(iResult).NumberOfSamples <= 2));
                
                % === OPTICAL FLOW ===
                if ~isempty(GlobalData.DataSet(iDS).Results(iResult).OpticalFlow)
                    sSurf = bst_memory('GetSurface', TessInfo(iTess).SurfaceFile);
                    panel_opticalflow('PlotOpticalFlow', hFig, GlobalData.DataSet(iDS).Results(iResult).OpticalFlow, ...
                                      GlobalData.UserTimeWindow.CurrentTime, sSurf); 
                end

            case 'Timefreq'
                % === LOAD TIMEFRQ VALUES ===
                % Get results index
                iTimefreq = bst_memory('GetTimefreqInDataSet', iDS, TessInfo(iTess).DataSource.FileName);
                % If Results file is not found in GlobalData structure
                if isempty(iTimefreq)
                    isOk = 0;
                    return
                end
                
                % === GET CURRENT VALUES ===
                % Get figure values
                TfInfo = getappdata(hFig, 'Timefreq');
                % Get results values
                TessInfo(iTess).Data = bst_memory('GetTimefreqValues', iDS, iTimefreq, TfInfo.RowName, TfInfo.iFreqs, 'CurrentTimeIndex', TfInfo.Function);
                % If min/max values for this file were not computed yet
                if isempty(TessInfo(iTess).DataMinMax)
                    TessInfo(iTess).DataMinMax = bst_memory('GetTimefreqMaximum', iDS, iTimefreq, TfInfo.Function);
                end
                % Reset Overlay cube
                TessInfo(iTess).OverlayCube = [];

                % Check the consistency between the number of results points (number of sources)
                % and the number of vertices of the target surface patch (IGNORE TEST FOR MRI)
                if strcmpi(TessInfo(iTess).Name, 'Anatomy')
                    % Nothing to check right now
                elseif ~isempty(TessInfo(iTess).DataSource.Atlas) 
                    if (length(TessInfo(iTess).Data) ~= length(TessInfo(iTess).DataSource.Atlas.Scouts))
                        bst_error(sprintf(['Number of sources (%d) is different from number of scouts (%d).\n\n' ...
                                  'Please compute the sources again.'], size(TessInfo(iTess).Data, 1), length(TessInfo(iTess).DataSource.Atlas.Scouts)), 'Data mismatch', 0);
                        isOk = 0;
                        return;
                    end
                elseif (length(TessInfo(iTess).Data) ~= TessInfo(iTess).nVertices) && ~strcmpi(TessInfo(iTess).Name, 'Anatomy')
                    bst_error(sprintf(['Number of sources (%d) is different from number of vertices (%d).\n\n' ...
                              'Please compute the sources again.'], size(TessInfo(iTess).Data, 1), TessInfo(iTess).nVertices), 'Data mismatch', 0);
                    isOk = 0;
                    return;
                end
                % Update "Static" status for this figure
                setappdata(hFig, 'isStatic',     (GlobalData.DataSet(iDS).Timefreq(iTimefreq).NumberOfSamples <= 2));
                setappdata(hFig, 'isStaticFreq', (size(GlobalData.DataSet(iDS).Timefreq(iTimefreq).TF,3) <= 1));
                
            case 'Dipoles'
                % === LOAD DIPOLES VALUES ===
                % Get results index
                iDipoles = bst_memory('GetDipolesInDataSet', iDS, TessInfo(iTess).DataSource.FileName);
                % If Results file is not found in GlobalData structure
                if isempty(iDipoles)
                    % Load Results file
                    [iDS, iDipoles] = bst_memory('LoadDipolesFile', TessInfo(iTess).DataSource.FileName);
                    if isempty(iDipoles)
                        return
                    end
                end
                
                % === GET CURRENT VALUES ===
                % Get results values
                % TessInfo(iTess).Data = bst_memory('GetDipolesValues', iDS, iDipoles, 'CurrentTimeIndex');
                % If min/max values for this file were not computed yet
                if isempty(TessInfo(iTess).DataMinMax)
                    TessInfo(iTess).DataMinMax = [0 100];
                end
                % Reset Overlay cube
                TessInfo(iTess).OverlayCube = [];

                % Update "Static" status for this figure
                setappdata(hFig, 'isStatic', (GlobalData.DataSet(iDS).Dipoles(iDipoles).NumberOfSamples <= 2));
                
            case 'Surface'
                % Get loaded surface
                SurfaceFile = TessInfo(iTess).DataSource.FileName;
                sSurf = bst_memory('LoadSurface', SurfaceFile);
                % Build uniform data vector
                TessInfo(iTess).Data = ones(length(sSurf.Vertices),1);
                TessInfo(iTess).DataMinMax = [.5 .5];
                setappdata(hFig, 'isStatic', 1);
                
            otherwise
                % Nothing to do
        end
        % Error if all data values are null
        if (max(abs(TessInfo(iTess).DataMinMax)) == 0)
            disp('BST> All values are null. Please check your input file.');
        end
    end
    % Update surface definition
    setappdata(hFig, 'Surface', TessInfo);
    % Update colormap
    UpdateSurfaceColormap(hFig, iSurfaces);
end



%% ===== UPDATE SURFACE COLORMAP =====
function UpdateSurfaceColormap(hFig, iSurfaces)
% disp('=== panel_surface > UpdateSurfaceColormap ===');
    % Get surfaces list 
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo)
        return
    end
    % If the aim is to update all the surfaces 
    if (nargin < 2) || isempty(iSurfaces)
        iSurfaces = 1:length(TessInfo);
    end
    
    % Get default colormap to use for this figure
    ColormapInfo = getappdata(hFig, 'Colormap');
    % Get figure axes
    hAxes = [findobj(hFig, '-depth', 1, 'Tag', 'Axes3D'), ...
             findobj(hFig, '-depth', 1, 'Tag', 'axc'), ...
             findobj(hFig, '-depth', 1, 'Tag', 'axa'), ...
             findobj(hFig, '-depth', 1, 'Tag', 'axs')];
    
    % Get figure index (in DataSet structure)
    [tmp__, iFig, iDS] = bst_figures('GetFigure', hFig);
    % Find the DataSet indice that corresponds to the current figure
    if isempty(iDS)
        error('No DataSet acessible for this 3D figure');
    end
    DataType = [];
    
    % For each surface
    for iTess = iSurfaces
        % ===== COLORMAPPING =====
        % Get colormap
        sColormap = bst_colormaps('GetColormap', TessInfo(iTess).ColormapType);

        % === Colormap : Normalized or Absolute ?
        % Dipoles density: only percentage
        if strcmpi(TessInfo(iTess).DataSource.Type, 'Dipoles')
            absMaxVal = 100;
        % Normalized : Color bounds (CLim) are set to a local extrema (at this time frame)
        elseif sColormap.isNormalized
            absMaxVal = max(abs(TessInfo(iTess).Data(:)));
        % Fixed : Color bounds are fixed, defined statically by the user (MaxValue)
        elseif ~isempty(sColormap.MaxValue)
            absMaxVal = sColormap.MaxValue;
        % Absolute : Color bounds (CLim) are set to the global extrema (over all the time frames)
        else
            absMaxVal = max(abs(TessInfo(iTess).DataMinMax));
        end
            
        % === Data values : Absolute | Normal ? ===
        if ismember(TessInfo(iTess).DataSource.Type, {'Dipoles'})
            TessInfo(iTess).DataLimitValue = [0, absMaxVal];
        elseif sColormap.isAbsoluteValues
            % Display absolute values of data
            TessInfo(iTess).Data = abs(TessInfo(iTess).Data);
            TessInfo(iTess).DataLimitValue = [0, absMaxVal];
        else
            % Display normal data (positive and negative values)
            TessInfo(iTess).DataLimitValue = [-absMaxVal, absMaxVal];
        end
        
        % If current colormap is the default colormap for this figure (for colorbar)
        if strcmpi(ColormapInfo.Type, TessInfo(iTess).ColormapType) && ~isempty(TessInfo(iTess).DataSource.FileName)
            if all(~isnan(TessInfo(iTess).DataLimitValue)) && (TessInfo(iTess).DataLimitValue(1) < TessInfo(iTess).DataLimitValue(2))
                set(hAxes, 'CLim', TessInfo(iTess).DataLimitValue);
            else
                %warning('Brainstorm:AxesError', 'Error using set: Bad value for axes property: CLim: Values must be increasing and non-NaN.');
                set(hAxes, 'CLim', [0 1e-30]);
            end
            DataType = TessInfo(iTess).DataSource.Type;
            % sLORETA: Do not use regular source scaling (pAm)
            isSLORETA = strcmpi(DataType, 'Source') && ~isempty(strfind(lower(TessInfo(iTess).DataSource.FileName), 'sloreta'));
            if isSLORETA 
                DataType = 'sLORETA';
            end
        end
                            
        % ===== DISPLAY ON MRI =====
        if strcmpi(TessInfo(iTess).Name, 'Anatomy') && ~isempty(TessInfo(iTess).DataSource.Type) && isempty(TessInfo(iTess).OverlayCube)
            % Progress bar
            isProgressBar = bst_progress('isVisible');
            bst_progress('start', 'Display MRI', 'Updating values...');
            % ===== Update surface display =====  
            % Update figure's appdata (surface list)
            setappdata(hFig, 'Surface', TessInfo);
            % Update OverlayCube
            TessInfo = UpdateOverlayCube(hFig, iTess);
            % Hide progress bar
            if ~isProgressBar
                bst_progress('stop');
            end
            % Put focus back on previous figure
            curFig = bst_figures('GetCurrentFigure', '3D');
            if ~isempty(curFig)
                figure(curFig);
            end
        else
            % Update figure's appdata (surface list)
            setappdata(hFig, 'Surface', TessInfo);
            % Update surface color
            figure_callback(hFig, 'UpdateSurfaceColor', hFig, iTess);
        end
    end
    
    % ===== Colorbar ticks and labels =====
    if ~isempty(ColormapInfo.Type)
        % Get figure colormap
        sColormap = bst_colormaps('GetColormap', ColormapInfo.Type);
        % Set figure colormap
        set(hFig, 'Colormap', sColormap.CMap);
        % Create/Delete colorbar
        bst_colormaps('SetColorbarVisible', hFig, sColormap.DisplayColorbar);
        % Display only one colorbar (preferentially the results colorbar)
        bst_colormaps('ConfigureColorbar', hFig, ColormapInfo.Type, DataType);
    % No colorbar should be displayed in this figure
    else
        % Delete colorbar
        bst_colormaps('SetColorbarVisible', hFig, 0);
    end
end

    
%% ===== GET SELECTED SURFACE =====
% Usage:  [iTess, TessInfo, sSurf, hFig] = GetSelectedSurface()
%         [iTess, TessInfo, sSurf, hFig] = GetSelectedSurface(hFig)
function [iTess, TessInfo, hFig, sSurf] = GetSelectedSurface(hFig) %#ok<DEFNU>
    % If target figure is not defined: use the current 3D figure
    if ((nargin < 1) || isempty(hFig))
        % Get current 3d figure
        hFig = bst_figures('GetCurrentFigure', '3D');
        % No current 3D figure: error
        if isempty(hFig)
            return
        end
    end
    % Get surface descriptions
    TessInfo = getappdata(hFig, 'Surface');
    iTess    = getappdata(hFig, 'iSurface');
    % Get the loaded structure in memory
    sSurf = [];
    if (nargout >= 4) && ~isempty(TessInfo) && ~isempty(iTess)
        sSurf = bst_memory('GetSurface', TessInfo(iTess).SurfaceFile);
    end
end


%% ===== GET SURFACE: ANATOMY =====
function [sMri,TessInfo,iTess,iMri] = GetSurfaceMri(hFig)
    sMri  = [];
    iTess = [];
    iMri  = [];
	% Get list of surfaces for the figure
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo)
        return
    end
    % Find "Anatomy"
    iTess = find(strcmpi({TessInfo.Name}, 'Anatomy'));
    if isempty(iTess)
        return
    elseif (length(iTess) > 1)
        iTess = iTess(1);
    end
    % Get Mri filename
    MriFile = TessInfo(iTess).SurfaceFile;
    % Get loaded MRI
    [sMri,iMri] = bst_memory('GetMri', MriFile);
end


%% ===== REMOVE A SURFACE =====
function RemoveSurface(hFig, iSurface)
    % Get figure appdata (surfaces configuration)
    TessInfo = getappdata(hFig, 'Surface');
    if (iSurface < 0) || (iSurface > length(TessInfo))
        return;
    end
    % Remove associated patch
    iRemPatch = ishandle(TessInfo(iSurface).hPatch);
    delete(TessInfo(iSurface).hPatch(iRemPatch));
    % Remove surface from the figure's surfaces list
    TessInfo(iSurface) = [];
    % Update figure's surfaces list
    setappdata(hFig, 'Surface', TessInfo);
    % Set another figure as current figure
    if isempty(TessInfo)
        setappdata(hFig, 'iSurface', []);
    elseif (iSurface <= length(TessInfo))
        setappdata(hFig, 'iSurface', iSurface);
    else
        setappdata(hFig, 'iSurface', iSurface - 1);
    end
    % Reload scouts
    panel_scouts('ReloadScouts', hFig);
end
       


%% ===== PLOT MRI =====
% Usage:  hs = panel_surface('PlotMri', hFig, posXYZ) : Set the position of cuts and plot MRI
%         hs = panel_surface('PlotMri', hFig)         : Plot MRI for current positions
function hs = PlotMri(hFig, posXYZ)
    % Get MRI
    [sMri,TessInfo,iTess,iMri] = GetSurfaceMri(hFig);
    % Set positions or use default
    if (nargin < 2) || isempty(posXYZ)
        posXYZ = TessInfo(iTess).CutsPosition;
        iDimPlot = ~isnan(posXYZ);
    else
        iDimPlot = ~isnan(posXYZ);
        TessInfo(iTess).CutsPosition(iDimPlot) = posXYZ(iDimPlot);
    end
    % Get initial threshold value
    threshold = TessInfo(iTess).SurfSmoothValue * 2 * double(sMri.Histogram.bgLevel);
    % Get colormaps
    sColormapData = bst_colormaps('GetColormap', TessInfo(iTess).ColormapType);
    sColormapMri  = bst_colormaps('GetColormap', 'anatomy');
    MriOptions = bst_get('MriOptions');
    % Define OPTIONS structure
    OPTIONS.sMri             = sMri;
    OPTIONS.iMri             = iMri;
    OPTIONS.cutsCoords       = posXYZ;                         % [x,y,z] location of the cuts in the volume
    OPTIONS.MriThreshold     = threshold;                      % MRI threshold (if value<threshold : background)
    OPTIONS.MriAlpha         = TessInfo(iTess).SurfAlpha;      % MRI alpha value (ie. opacity)
    OPTIONS.MriColormap      = sColormapMri.CMap;              % MRI Colormap     
    OPTIONS.OverlayCube      = TessInfo(iTess).OverlayCube;    % Overlay values
    OPTIONS.OverlayThreshold = TessInfo(iTess).DataThreshold;  % Overlay threshold
    OPTIONS.OverlaySizeThreshold = TessInfo(iTess).SizeThreshold;  % Overlay size threshold
    OPTIONS.OverlayAlpha     = TessInfo(iTess).DataAlpha;      % Overlay transparency
    OPTIONS.OverlayColormap  = sColormapData.CMap;             % Overlay colormap
    OPTIONS.OverlayBounds    = TessInfo(iTess).DataLimitValue; % Overlay colormap amplitude, [minValue,maxValue]
    OPTIONS.isMipAnatomy     = MriOptions.isMipAnatomy;
    OPTIONS.isMipFunctional  = MriOptions.isMipFunctional;
    OPTIONS.MipAnatomy       = TessInfo(iTess).MipAnatomy;
    OPTIONS.MipFunctional    = TessInfo(iTess).MipFunctional;
    % Plot cuts
    [hs, OutputOptions] = mri_draw_cuts(hFig, OPTIONS);

    TessInfo(iTess).hPatch(iDimPlot) = hs(iDimPlot);
    % Save maximum in each direction in TessInfo structure
    if OPTIONS.isMipAnatomy
        iUpdateSlice = ~cellfun(@isempty, OutputOptions.MipAnatomy);
        TessInfo(iTess).MipAnatomy(iUpdateSlice) = OutputOptions.MipAnatomy(iUpdateSlice);
    end
    if OPTIONS.isMipFunctional
        iUpdateSlice = ~cellfun(@isempty, OutputOptions.MipFunctional);
        TessInfo(iTess).MipFunctional(iUpdateSlice) = OutputOptions.MipFunctional(iUpdateSlice);
    end
    % Save TessInfo
    setappdata(hFig, 'Surface', TessInfo);
end


%% ===== UPDATE OVERLAY MASKS =====
function UpdateOverlayCubes(hFig) %#ok<DEFNU>
    for i = 1:length(hFig)
        [sMri, TessInfo, iTess] = GetSurfaceMri(hFig(i));
        if ~isempty(iTess) && ~isempty(TessInfo(iTess).Data)
            UpdateOverlayCube(hFig(i), iTess);
        end
    end
end


%% ===== UPDATE OVERLAY MASK =====
% Usage:  TessInfo = UpdateOverlayCube(hFig, iTess)
function TessInfo = UpdateOverlayCube(hFig, iTess)
    global GlobalData;
    % Get MRI
    TessInfo = getappdata(hFig, 'Surface');
    sMri = bst_memory('GetMri', TessInfo(iTess).SurfaceFile);
    if isempty(sMri) || isempty(sMri.Cube) || (isempty(TessInfo(iTess).Data) && ~strcmpi(TessInfo(iTess).DataSource.Type, 'Dipoles'))
       return 
    end
    isProgressBar = 0;
    SurfaceFile = [];
    OverlayCube = [];
    isVolumeGrid = 0;
    % Process depend on overlay data file
    switch (TessInfo(iTess).DataSource.Type)
       case 'Data'
           % Get scalp surface
           error('Not supported yet');
        case 'Source'
            % Get cortex surface
            sSubject = bst_get('MriFile', sMri.FileName);
            SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
            % Get loaded results file
            [iDS, iResult] = bst_memory('GetDataSetResult', TessInfo(iTess).DataSource.FileName);
            if isempty(iDS)
                return
            end            
            % Check source grid type
            isVolumeGrid = strcmpi(GlobalData.DataSet(iDS).Results(iResult).HeadModelType, 'volume');
        case 'Timefreq'
            % Get cortex surface
            sSubject = bst_get('MriFile', sMri.FileName);
            SurfaceFile = sSubject.Surface(sSubject.iCortex).FileName;
            % Get loaded timefreq file
            [iDS, iTf] = bst_memory('GetDataSetTimefreq', TessInfo(iTess).DataSource.FileName);
            if isempty(iDS)
                return
            end
            % If timefreq on sources
            if strcmpi(GlobalData.DataSet(iDS).Timefreq(iTf).DataType, 'results')
                % Check source grid type
                if ~isempty(GlobalData.DataSet(iDS).Timefreq(iTf).DataFile)
                    [iDS, iResult] = bst_memory('GetDataSetResult', GlobalData.DataSet(iDS).Timefreq(iTf).DataFile);
                    isVolumeGrid = strcmpi(GlobalData.DataSet(iDS).Results(iResult).HeadModelType, 'volume');
                else
                    isVolumeGrid = 0;
                end
            end
        case 'Surface'
            % Get surface specified in DataSource.FileName
            SurfaceFile = TessInfo(iTess).DataSource.FileName;
        case 'Dipoles'
            sDipoles = panel_dipoles('GetSelectedDipoles', hFig);
            OverlayCube = panel_dipoles('ComputeDensity', sMri, sDipoles);
    end
    
    
    % ===== 3D OVERLAYS =====
    if ~isempty(OverlayCube)
        TessInfo(iTess).OverlayCube = OverlayCube;
    % ===== DISPLAY SURFACE/GRIDS =====
    else
        % === INTERPOLATION MRI<->SURFACE ===
        [sSurf, iSurf] = bst_memory('LoadSurface', SurfaceFile);
        tess2mri_interp = bst_memory('GetTess2MriInterp', iSurf);
        % If no interpolation tess<->mri accessible : exit
        if isempty(tess2mri_interp)
           return 
        end
        % Progress bar
        isProgressBar = bst_progress('isVisible');
        bst_progress('start', 'Display MRI', 'Updating values...');
        % === INTERPOLATION MRI<->GRID ===
        if isVolumeGrid
            % Compute interpolation
            MriInterp = bst_memory('GetGrid2MriInterp', iDS, iResult);
        else
            % Only surface interpolation is needed
            MriInterp = tess2mri_interp;
        end
        % === ATLAS SOURCES ===
        if ~isempty(TessInfo(iTess).DataSource.Atlas)
            % Initialize full cortical map
            DataScout = TessInfo(iTess).Data;
            DataSurf = zeros(size(MriInterp,2),1);
            % Duplicate the value of each scout to all the vertices
            sScouts = TessInfo(iTess).DataSource.Atlas.Scouts;
            for i = 1:length(sScouts)
                DataSurf(sScouts(i).Vertices,:) = DataScout(i,:);
            end
        else
            DataSurf = TessInfo(iTess).Data;
        end
        % === UPDATE MASK ===
        mriSize = size(sMri.Cube);
        % Build interpolated cube
        TessInfo(iTess).OverlayCube = tess_interp_mri_data(MriInterp, mriSize, DataSurf, isVolumeGrid);
    end
    
    % === UPDATE DISPLAY ===
    % Reset MIP functional fields
    TessInfo(iTess).MipFunctional = cell(3,1);    
    % Get surface description
    setappdata(hFig, 'Surface', TessInfo);
    % Redraw surface vertices color
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iTess);
    % Hide progress bar
    if ~isProgressBar
        bst_progress('stop');
    end
end


%% ===== SET SURFACE TRANSPARENCY =====
function SetSurfaceTransparency(hFig, iSurf, alpha)
    % Update surface transparency
    TessInfo = getappdata(hFig, 'Surface');
    TessInfo(iSurf).SurfAlpha = alpha;
    setappdata(hFig, 'Surface', TessInfo);
    % Update panel controls
    UpdateSurfaceProperties();
    % Update surface display
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurf);
    figure_callback(hFig, 'UpdateSurfaceAlpha', hFig, iSurf);
end


%% ===== SET DATA THRESHOLD =====
function SetDataThreshold(hFig, iSurf, value) %#ok<DEFNU>
    % Update surface transparency
    TessInfo = getappdata(hFig, 'Surface');
    TessInfo(iSurf).DataThreshold = value;
    setappdata(hFig, 'Surface', TessInfo);
    % Update panel controls
    UpdateSurfaceProperties();
    % Update color display on the surface
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurf);
end

%% ===== SET SIZE THRESHOLD =====
function SetSizeThreshold(hFig, iSurf, value) %#ok<DEFNU>
    % Update surface transparency
    TessInfo = getappdata(hFig, 'Surface');
    TessInfo(iSurf).SizeThreshold = value;
    setappdata(hFig, 'Surface', TessInfo);
    % Update panel controls
    UpdateSurfaceProperties();
    % Update color display on the surface
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurf);
end

%% ===== SET SURFACE SMOOTH =====
function SetSurfaceSmooth(hFig, iSurf, value)
    % Update surface transparency
    TessInfo = getappdata(hFig, 'Surface');
    TessInfo(iSurf).SurfSmoothValue = value;
    setappdata(hFig, 'Surface', TessInfo);
    % Update panel controls
    UpdateSurfaceProperties();
    % For MRI display : Smooth slider changes threshold
    if strcmpi(TessInfo(iSurf).Name, 'Anatomy')
        figure_callback(hFig, 'UpdateMriDisplay', hFig, [], TessInfo, iSurf);
    % Else: Update color display on the surface
    else
        % Smooth surface
        figure_callback(hFig, 'UpdateSurfaceAlpha', hFig, iSurf);
        % Update scouts displayed on this surfce
        panel_scouts('UpdateScoutsVertices', TessInfo(iSurf).SurfaceFile);
        % Set the new value as the default value
        DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
        DefaultSurfaceDisplay.SurfSmoothValue = TessInfo(iSurf).SurfSmoothValue;
        bst_set('DefaultSurfaceDisplay', DefaultSurfaceDisplay);
    end
end


%% ===== SET AUTO TRANSPARENCY =====
function SetAutoTransparency(hFig)
    % Get surfaces definitions
    TessInfo = getappdata(hFig, 'Surface');
    % Look for different surfaces types
    iOrder = [find(ismember({TessInfo.Name}, {'Cortex','Anatomy'})), ...
              find(strcmpi({TessInfo.Name}, 'InnerSkull')), ...
              find(strcmpi({TessInfo.Name}, 'OuterSkull')), ...
              find(strcmpi({TessInfo.Name}, 'Scalp'))];
    % Set other surfaces transparency if cortex at the same time
    for i = 2:length(iOrder)
        SetSurfaceTransparency(hFig, iOrder(i), 0.7);
    end
end
    

%% ===== SET SURFACE COLOR =====
function SetSurfaceColor(hFig, iSurf, colorCortex, colorSulci)
    % Compute the color used to display sulci
    if (nargin < 4) || isempty(colorSulci)
        colorSulci = .73 .* colorCortex;
    end
    % Get description of surfaces
    TessInfo = getappdata(hFig, 'Surface');
    % Update surface description (figure's appdata)
    TessInfo(iSurf).AnatomyColor(1,:) = colorSulci;
    TessInfo(iSurf).AnatomyColor(2,:) = colorCortex;
    % Update Surface appdata structure
    setappdata(hFig, 'Surface', TessInfo);
    
    % Get panel controls
    ctrl = bst_get('PanelControls', 'Surface');
    % Change button color
    ctrl.jButtonSurfColor.setBackground(java.awt.Color(colorCortex(1), colorCortex(2), colorCortex(3)));
    % Update panel controls
    UpdateSurfaceProperties();
    
    % Update color display on the surface
    figure_callback(hFig, 'UpdateSurfaceColor', hFig, iSurf);
end

%% ===== APPLY DEFAULT DISPLAY TO SURFACE =====
function ApplyDefaultDisplay() %#ok<DEFNU>
    % Get panel controls
    ctrl = bst_get('PanelControls', 'Surface');
    % Get defaults for surface display
    DefaultSurfaceDisplay = bst_get('DefaultSurfaceDisplay');
    % Surface smooth
    if (ctrl.jSliderSurfSmoothValue.getValue() ~= DefaultSurfaceDisplay.SurfSmoothValue * 100)
        ctrl.jSliderSurfSmoothValue.setValue(DefaultSurfaceDisplay.SurfSmoothValue * 100);
        event = java.awt.event.MouseEvent(ctrl.jSliderSurfSmoothValue, 0, 0, 0, 0, 0, 1, 0, 0);
        SliderCallback([], event, 'SurfSmoothValue');
    end
    % Surface edges
    if DefaultSurfaceDisplay.SurfShowSulci && ~ctrl.jButtonSurfSulci.isSelected()
        ctrl.jButtonSurfSulci.doClick();
    end
end

