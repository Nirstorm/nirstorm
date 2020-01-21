function [varargout] = be_viz_result(RESULTS, varargin)

% ===== Process inputs =====
opacite     =   .5;
vue         =   [0 90];
chx         =   cellfun( @(a) num2str(a), mat2cell(0:.05:1,1,ones(1,21)), 'uni', 0 );
for ii   =   1 : numel(varargin)
    
    if isnumeric(varargin{ii})
        varargin{ii}    =   num2str(varargin{ii});
    end
    
    switch varargin{ii}
        
        case 'hemiL'
            SIDE    =   'left';
            
        case 'hemiR'
            SIDE    =   'right';
            
        case 'top'
            vue     =   [0 90];
            
        case 'bottom'
            vue     =   [-90 90];
            
        case 'right'
            vue     =   [0 0];
            
        case 'left'
            vue     =   [180 0];
            
        case chx
            opacite =   str2double(varargin{ii});
            fprintf('\nBE_VIZ_RESULT>>\t transparency set to %1.1f\n', opacite);
    end
    
end

% ===== CREATE NEW TEMPORARY FILE =====
iProt               =   bst_get('ProtocolInfo');
[sStudy,iStudy]     =   bst_get('Study');
[sSubject]          =   bst_get('Subject', sStudy.BrainStormSubject);
SURFACE             =   sSubject.Surface(sSubject.iCortex).FileName;
[a,b,c]             =   bst_fileparts( sStudy.FileName );
FLNM                =   be_fullfile(a,'results_TEMPORARY.mat');
TEMP                =   struct;
TEMP.ImageGridAmp   =   [RESULTS(:,1) RESULTS];
TEMP.ImagingKernel  =   [];
TEMP.Time           =   [0 0:(size(RESULTS,2)-1)];
TEMP.nComponents    =   1;
TEMP.SurfaceFile    =   SURFACE;
TEMP.Comment        =   'This file is temporary and must be deleted';
TEMP.HeadModelType  =   'surface';
TEMP.DataFile       =   '';
save( be_fullfile(iProt.STUDIES,FLNM), '-struct', 'TEMP');

% ===== SET STUDY FILE =====
newResult               = db_template('Results');
newResult.Comment       = 'This file is temporary and must be deleted';
newResult.FileName      = FLNM;
newResult.DataFile      = '';
newResult.isLink        = 0;
newResult.HeadModelType = 'Surface';
iResult                 = length(sStudy.Result) + 1;
sStudy.Result(iResult)  = newResult;
bst_set('Study', iStudy, sStudy);
db_links('Study', iStudy);
% panel_protocols('UpdateNode', 'Study', iStudy);
% db_update(999);
% db_reload_studies(iStudy);

% ===== SET GLOBAL VARIABLE
global GlobalData
nDS     =   numel(GlobalData.DataSet)+1;
GlobalData.DataSet(nDS).DataFile    =   '';
GlobalData.DataSet(nDS).StudyFile   =   sStudy.FileName;
GlobalData.DataSet(nDS).SubjectFile =   sStudy.BrainStormSubject;
GlobalData.DataSet(nDS).ChannelFile =   sStudy.Channel.FileName;
GlobalData.DataSet(nDS).Surfaces    =   struct;
GlobalData.DataSet(nDS).Measures    =   struct('DataType', [], 'F', [], 'Time', TEMP.Time, 'SamplingRate', 0, 'NumberOfSamples', size(TEMP.Time,2), 'ChannelFlag', [], 'sFile', [], 'isModified', 0, 'ColormapType', '');
GlobalData.DataSet(nDS).Results     =   TEMP;
GlobalData.DataSet(nDS).Results.FileName    =   FLNM;
GlobalData.DataSet(nDS).Results.DataType    =   'results';
GlobalData.DataSet(nDS).Results.SamplingRate    =   0;
GlobalData.DataSet(nDS).Results.NumberOfSamples =   size(TEMP.Time,2);
GlobalData.DataSet(nDS).Results.ChannelFlag     =   [];
GlobalData.DataSet(nDS).Results.ColormapType    =   '';
GlobalData.DataSet(nDS).Results.Atlas           =   '';
GlobalData.DataSet(nDS).Results.OpticalFlow     =   '';
db_update(999);

% ===== SET UP NEW FIGURE =====
FigureId = db_template('FigureId');
FigureId.Type     = '3DViz';
FigureId.SubType  = '';
FigureId.Modality = [];

% Create figure
[hFig] = bst_figures('CreateFigure', nDS, FigureId, 'AlwaysCreate');
iSurf = panel_surface('AddSurface', hFig, SURFACE);
panel_surface('SetSurfaceSmooth', hFig, iSurf, .5);
setappdata(hFig, 'StudyFile',    sStudy.FileName);
setappdata(hFig, 'DataFile',     FLNM);
setappdata(hFig, 'SubjectFile',  sStudy.BrainStormSubject);
TessInfo    =   getappdata(hFig, 'Surface');
TessInfo(iSurf).SurfAlpha = opacite;
if exist('SIDE', 'var')
    TessInfo(iSurf).Resect = SIDE;
end
setappdata(hFig, 'Surface', TessInfo);
figure_3d('UpdateSurfaceAlpha', hFig, iSurf);


% ===== OVERLAY DATA ON SURFACE =====
% Set data source for this surface
isOk = panel_surface('SetSurfaceData', hFig, iSurf, 'Source', FLNM, 0);
panel_surface('ApplyDefaultDisplay');
figure_3d('SetStandardView', hFig, 'top');
set(hFig,'Visible', 'on');
delete( be_fullfile(iProt.STUDIES,FLNM) )
sStudy.Result(iResult)  = [];
bst_set('Study', iStudy, sStudy);
db_links('Study', iStudy);
GlobalData.DataSet(nDS) = [];
panel_surface('SetDataThreshold', hFig, iSurf, 0);
set(gca, 'view', vue);

if nargout>0
    varargout{1}    =   hFig;
end
if nargout>1
    varargout{2}    =   findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
end

end