function [hFig] = be_viz_topo(channels, data)

global GlobalData


% ===== SET UP NEW FIGURE =====
FigureId            =   db_template('FigureId');
FigureId.Type       =   '3DViz';
FigureId.SubType    =   '';
FigureId.Modality   =   [];

% Create figure
hFig                = figure_3d('CreateFigure', FigureId);
hAxes               =   findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
hChildren           = get(hAxes, 'Children');
	delete(hChildren(~strcmpi(get(hChildren, 'Type'), 'light')));
	% Set Topography axes as current axes
	set(0,    'CurrentFigure', hFig);
	set(hFig, 'CurrentAxes',   hAxes);
	hold on
	% Reset initial zoom factor and camera position
	figure_3d('ResetView', hFig);
    
    
% ===== MAKE TOPOGRAPHY =====
% commons
chan_loc        =   ( cellfun( @(a) mean(a(:,1:4),2)', {channels.Loc}, 'uni', 0) )';
chan_loc        =   cell2mat( chan_loc );

% Compute best fitting sphere from sensors
[bfs_center, bfs_radius] = bst_bfs(chan_loc);


% ===== CREATE A HIGH-DEF SURFACE =====
% Remove the duplicated positions
precision = 1e5;
Vertices = unique(round(chan_loc * precision)/precision,'rows');

% Tesselate sensor cap
Faces = channel_tesselate(Vertices, 1);

% Clean from some very pathological triangles
Faces = tess_threshold(Vertices, Faces, 3, []);

% Refine mesh
[Vertices, Faces] = tess_refine(Vertices, Faces, [], 1);
if (length(Vertices) < 800)
    [Vertices, Faces] = tess_refine(Vertices, Faces, [], 1);
end


% ===== MAKE INTERPOLATION =====
% Get extrapolation options
epsilon                 =   1e-4;
Whitener                =   [];

% Compute normals at each vertex
VerticesNormals         =   tess_normals(Vertices, Faces);

% Create a pseudo-channel file of virtual magnetometers
destChan                =   repmat(db_template('channeldesc'), [1, length(Vertices)]);
for i = 1:length(Vertices)
    destChan(i).Name    =   num2str(i);
    destChan(i).Type    =   'MEG';
    destChan(i).Loc     =   Vertices(i,:)';
    destChan(i).Orient  =   VerticesNormals(i,:)';
    destChan(i).Weight  =   1;
end

% Compute interpolation
[WExtrap, src_xyz]      = 	be_computeinterp(channels, destChan, bfs_center, bfs_radius, Whitener, epsilon);

% Apply interpolation matrix sensors => display surface
DataToPlot              =   WExtrap * data;


% ===== DRAW INTERPOLATION =====
% Get figure colormap
ColormapInfo            =   getappdata(hFig, 'Colormap');
sColormap               =   bst_colormaps('GetColormap', ColormapInfo.Type);
    
% ===== Map data on target patch =====
hSurf = patch('Faces',     Faces, ...
    'Vertices',  Vertices, ...
    'BackfaceLighting', 'lit', ...
    'Parent',    hAxes, ...
    'FaceVertexCData', DataToPlot, ...
    'EdgeColor', 'none', ...
    'FaceColor', 'interp');
           
% Surface lighting
material ([.95 0 0])
lighting gouraud
set(hFig, 'Visible', 'on')

return