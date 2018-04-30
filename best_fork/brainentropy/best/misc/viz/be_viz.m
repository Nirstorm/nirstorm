function viz(varargin)
% This function displays the continuous time-frequency plane and linked 
% ridges extracted using the RIDGEFILTER process. This function is  only
% operational under brainstorm environment. A study must be selected before
% launching VIZ.m
%
% Input:
% ------
%
%   []          :	Displays both time-frequency and ridges planes and 
%                   highlights ridge lines
%   w           :	Displays only time-frequency plane
%   r           	Display only ridges plane
%   zoom        :	Limits displayed planes to time window and frequency
%                   range of ridge analysis
%   -           :	Do not highlight ridge lines
%   'integer'   :   (vector) Highlight only ridges # 'integer'
%
%
% No output
% ---------
%
% -------------------------------------------------------------------------
%
% LATIS team, 2012
%
%% ==============================================
% License 
%
% BEst is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BEst is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BEst. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------


% Basic BST variables
iP  = bst_get('ProtocolInfo');


% Keeps variables in memory - faster execution
global vizData
if isfield(vizData, 'hpar')
    hpar = vizData.hpar;
    delete( get(hpar, 'Children') )
    TF = vizData.TF;
    TFfile = vizData.TF.TFfile;
else
    hpar = figure('DeleteFcn', 'clear global vizData');
    vizData.hpar = hpar;
    TFfile  = find_TF(iP);
    if isempty(TFfile)
        return
    end
    TF = load( be_fullfile(iP.STUDIES, TFfile), 'Time', 'Freqs', 'Options', 'ridgeMAP', 'DataTypes', 'TF' );
    TF.TFfile = TFfile;
    load( be_fullfile(iP.STUDIES, TF.TFfile), 'DataFile' )
    load( be_fullfile(iP.STUDIES, DataFile), 'Time')
    TF.rdgTIME = Time;
    TF.DataFile = DataFile;
    
    % Removes low freq ridges
    LIMs = be_closest(TF.Freqs, 1);
    for ii = 1:numel(TF.ridgeMAP{1})
        TF.ridgeMAP{1}{ii}(1:LIMs,:) = 0;
        TF.TF{1}{ii}(1:LIMs,:) = 0;
    end
    
end
vizData.lastcall = varargin;


%% GET RIDGE DATA
iD  = find_RDG(iP, TFfile);

% Process input arguments
dt    = diff(TF.Time([1 2]));
w = 0; r = 0; cx = 0; Tlims = fix( (TF.Time([1 end])-TF.Time(1)) / dt ) + 1; Flims = [1 numel(TF.Freqs)];
iDn = []; iDm = {}; rdgLIMS = [];
for ii = 1 : nargin
    switch varargin{ii}
        case '-'
            iD = [];
        case 'w'
            w = 1;
        case 'r'
            r = 1;
        case 'zoom'
            Tlims = TF.Options.RidgesWindow([1 end]);
            Tlims = fix( (Tlims - TF.Time(1)) / dt ) + 1;
            
            Flims    = TF.Options.RidgesFRange;
            [dum,Fl1]= min( abs(TF.Freqs-Flims(1)) );
            [dum,Fl2]= min( abs(TF.Freqs-Flims(2)) );
            Flims    = [Fl1 Fl2];
        case {'meg','eeg'}
            iDm     =   [iDm varargin{ii}]; 
        case 'Flims'
            try
                NMarg   =   eval(varargin{ii+1});
                if isnumeric(NMarg) && numel(NMarg)==2
                    Flims = be_closest(NMarg, TF.Freqs);
                end
            end
        case 'cx'
            cx = 1;
            try
                NMarg   =   eval(varargin{ii+1});
                if isnumeric(NMarg) && numel(NMarg)==2
                    rdgLIMS = NMarg;
                end
            end
        otherwise
            NMarg   =   str2double(varargin{ii});
            if ~isnan(NMarg) && isnumeric( NMarg ) && numel(NMarg)==1
                iDn = [iDn fix( str2double(varargin{ii}) )];
            end
    end
end
if ~(w+r)
    w = 1;
    r = 1;
end
if ~isempty(iDn)
    [dum, pnt] = ismember(iDn, [iD{:,2}] );
    iD = iD(pnt, :);
end
if isempty(iDm)
    iDm = TF.DataTypes;
end


% Show cortical map
iS	= bst_get('Study', be_fullfile(bst_fileparts(TF.DataFile),'brainstormstudy.mat') );
if cx == 1 && ~isempty( find_RSLT(iD, iS, iP, rdgLIMS) )
    iSS = bst_get('Subject', iS.BrainStormSubject);    
    TF.CXfile = iSS.Surface(iSS.iCortex);
    int = find_RSLT(iD, iS, iP, rdgLIMS);
    if ~isfield(vizData, 'CX')
        vizData.CX = load( be_fullfile(iP.SUBJECTS, TF.CXfile.FileName) );
    end
end

    
%% MAKE PLOTS
a_w = []; a_r = []; a_x = [];
for ii = 1 : numel(iDm)
    if w
        a_w(ii) = subplot(numel(iDm),w+r+cx,(ii-1)*(w+r+cx)+1, 'parent', hpar);
        show_TF( TF, 'cwt', a_w, ii )
        hold on
        adjust_window(TF, a_w, ii, Tlims, Flims, 1)
    end
    if r
        a_r(ii) = subplot(numel(iDm),w+r+cx,(ii-1)*(w+r+cx)+2, 'parent', hpar);
        show_TF( TF, 'rdg', a_r, ii )
        hold on
        adjust_window(TF, a_r, ii, Tlims, Flims, 0)
    end
    if cx
        a_x(ii) = subplot(numel(iDm),w+r+cx,(ii-1)*(w+r+cx)+3, 'parent', hpar);
        if ii==1 
        int = cellfun(@(a) mean(abs(a),2), int, 'uni', false);
        int = max([int{:}], [], 2);
        patch('Faces', vizData.CX.Faces, 'Vertices', vizData.CX.Vertices, 'EdgeColor', 'none', 'FaceColor', 'interp', 'FaceVertexCData', int);
        axis off square; daspect([1 1 1]); rotate3d on
        end
    end    
end

% Highlight ridges
show_R(iP, iD, rdgLIMS, [a_r a_w]);
    
   
return




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- HELPER FUNCTIONS --------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [TF] = find_TF(iP)

    % Get selected node
    jTreeProtocols  = bst_get('PanelElement', 'protocols', 'jTreeProtocols');
    selectedPaths   = awtinvoke(jTreeProtocols, 'getSelectionPaths()');
    lastnode        = awtinvoke(selectedPaths(1), 'getLastPathComponent()');

    % Process selected node
    FL  =   char( lastnode.getFileName() );
    TF  =   '';

    % Check node type - timefreq
    TP  =   strfind( be_fullfile(iP.STUDIES, FL), 'timefreq' );
    if numel(TP)==1
        O = load( be_fullfile(iP.STUDIES, FL), 'Options' );
        if strcmpi(O.Options.WaveletType, 'Morse')
            TF = FL;
        else
            fprintf('In viz.m:\tThis is not an analytic wavelet transform\n\tSelect a new TF file\n\n')
        end
        return
    end

    % Check node type - ridgedata
    TP  =   whos( '-file', be_fullfile(iP.STUDIES, FL), 'Fi' );
    if numel(TP)==1 
        O = load( be_fullfile(iP.STUDIES, FL), 'TFfile', 'DataType' );
        if exist( be_fullfile(iP.STUDIES, O.TFfile), 'file' ) && any(strcmpi(O.DataType, {'data', 'recordings'}))
            TF = O.TFfile;
        else
            fprintf('In viz.m:\tTF plane linked to this file was not found. Recompute\n\n')
        end
        return
    end
    
    % Check node type - data
    TP  =   whos( '-file', be_fullfile(iP.STUDIES, FL), 'DataType', 'F');
    if numel(TP)==2 && isempty(TF)
        iS       = bst_get('Study', be_fullfile(bst_fileparts(FL),'brainstormstudy.mat') );
        for ii = 1 : numel(iS.Timefreq)
            O = load( be_fullfile(iP.STUDIES, iS.Timefreq(ii).FileName), 'Options', 'DataFile' );
            if strcmpi(O.Options.WaveletType, 'Morse') && strcmpi( FL, O.DataFile )
                TF = iS.Timefreq(ii).FileName;
            end
        end
        return
    end
        
return
        
        
function [iD] = find_RDG(iP, TFfile)
    iS  =   bst_get( 'Study', be_fullfile( bst_fileparts(TFfile), 'brainstormstudy.mat' ) );
    iSd =   { iS.Data(:).FileName };
    iD  =   cell(0,2);
    
    warning('OFF')
    for ii = 1 : numel(iSd)
        a = load( be_fullfile(iP.STUDIES, iSd{ii}), 'TFfile' );
        if isfield(a, 'TFfile') && strcmp(TFfile, a.TFfile)
            rdgStr = strfind( iSd{ii}, 'RidgeFilter' );
            endPos = strfind( iSd{ii}, '.' );
            iD = [iD; iSd(ii) {str2double(iSd{ii}(rdgStr(end)+12:endPos-1))}];
        end
    end
            

return


function [int] = find_RSLT(iD, iS, iP, lm)
    int = {};
    dtF = {iS.Result.DataFile};
    rsF = {iS.Result.FileName};
    for ii =1 : size(iD,1)
        idx = find( strcmpi( iD{ii,1}, dtF ) );
        iii = load( be_fullfile( iP.STUDIES,rsF{idx} ), 'ImageGridAmp' );
        
        if ~isempty(lm)
            st = fix( size(iii.ImageGridAmp,2)*lm(1) ) + 1;
            nd = ceil( size(iii.ImageGridAmp,2)*lm(2) );
            int = [int {iii.ImageGridAmp(:,st:nd)}];
        else
            int = [int {iii.ImageGridAmp}];
        end     
    end    

return


function show_TF(TF, wtp, a, ii)

% Keep variables - faster execution
global vizData
if isfield(vizData, 'POW')
    POW{ii} = vizData.POW{ii};
    TF      = vizData.TF;
else
    POW{ii}             = abs(be_avgcell(TF.TF{ii})).^2;
    vizData.POW{ii}     = POW{ii};
    vizData.TF          = struct;
    vizData.TF.Time     = TF.Time;
    vizData.TF.Freqs    = TF.Freqs;
    vizData.TF.Options  = TF.Options;
    vizData.TF.ridgeMAP = TF.ridgeMAP;
    vizData.TF.TFfile   = TF.TFfile;
    vizData.TF.DataTypes= TF.DataTypes;
    vizData.TF.rdgTIME  = TF.rdgTIME;
    vizData.TF.DataFile = TF.DataFile;
end
    
% TF/ridge plane limits
TA    = TF.Time;
FA    = TF.Freqs;
Ytic  = 2.^(fix( log2(FA(1)) ):fix( log2(FA(end)) ));
temp  = num2str(Ytic','%4.2e');
Yticl = mat2cell( temp, ones(1, numel(Ytic)), size(temp,2) );

% Draw plane
switch wtp
    case 'cwt'
        imagesc(TA, log2(FA), POW{ii},'Parent',a(ii) );
        text(0, .5, TF.DataTypes{ii}, 'Parent',a(ii), 'Units', 'normalized');
    
    case 'rdg'
        TA    = TF.rdgTIME;
        imagesc(TA, log2(FA), be_avgcell(TF.ridgeMAP{ii}),'Parent',a(ii) );
end
set(a(ii), 'XLim', TA([1 end]), 'YGrid','on', 'YLim',log2(FA([1 end])),...
    'YTick',log2(Ytic(:)), 'YTickLabel',Yticl, 'YDir', 'normal', ...
    'FontName','Times', 'FontAngle','Italic', 'FontSize',8);
xlabel(a(ii),'time (s)','FontSize',8); 
ylabel(a(ii),'frequency (Hz)','FontSize',8);

return


function show_R(iP, iD, lm, a)

global vizData

% Loops on ridge lines
nbF    = size(vizData.POW{1}, 1);

for ii = 1 : size(iD,1)
    % get ridges points indices
    load( be_fullfile(iP.STUDIES, iD{ii,1}), 'TFpath' );
    x = fix(TFpath/nbF) + 1;
    y = mod(TFpath, nbF);
    if ~isempty(lm)
        st = fix( numel(TFpath)*lm(1) ) + 1;
        nd = ceil( numel(TFpath)*lm(2) );
        x = fix(TFpath(st:nd)/nbF) + 1;
        y = mod(TFpath(st:nd), nbF);
    end
    
    % Transform to plane coordinates
    xl = [0.95 1.05] .* vizData.TF.rdgTIME( [min(x) max(x)] );
    yl = log2( [0.95 1.05] .* vizData.TF.Freqs( [min(y) max(y)] ) );
    
    % Draw highlight boxes
    for jj = 1 : numel(a)
        axes(a(jj))
        FV.Vertices = [ xl(1) yl(1); xl(1) yl(2); xl(2) yl(2); xl(2) yl(1) ];
        FV.Faces    = [1 2 3 4 1];
        patch( FV, 'FaceColor', 'y', 'FaceAlpha', .5, 'LineWidth', 2)
        text( mean(xl), yl(2), sprintf('#%u (%.1fHz)', iD{ii,2}, median(vizData.TF.Freqs(y)) ), 'Fontsize', 20,...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','bottom' )
    end
end

return


function adjust_window(TF, a, ii, tlim, flim, col)

global vizData

% Get ridge analysis box coordinates
minC = min(min( ( vizData.POW{ii}( flim(1):flim(2), tlim(1) : tlim(end))) ) );
maxC = max(max( ( vizData.POW{ii}( flim(1):flim(2), tlim(1) : tlim(end))) ) );
set(a(ii), 'xlim', TF.Time(tlim));
set(a(ii), 'ylim', log2(TF.Freqs(flim)));
if col
    set(a(ii), 'clim', [minC maxC]);
end
    
return


function na = smoothabs(a)

[n,m,p] = size(a);
na = zeros(m,p);
for ii = 1 : n
    slice = squeeze( a(ii,:,:) );
    slice = abs( slice );
    na = na + slice;
end
    
return