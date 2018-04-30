function [NewFile] = be_store_me(W, varargin)

% What kind of data?
NewFile     =   'No file stored';
[kind,modality,iStudy,iAtlas]  =   what_kind(W);
if isempty( kind )
    fprintf('\nBE_STORE_ME error>> unrecognized data\n')
    return
end

% What parent?
parent  =   '';
if numel(varargin)>0 
    parent          =   varargin{1};   
end

% comment?
comment =   '';
if numel(varargin)>1 
    comment         =   varargin{2};   
end

% labels?
labels  =   '';
if numel(varargin)>2 
    labels          =   varargin{3};   
end

% found?
if isempty(kind)
    fprintf('\nBE_STORE_ME error>> unrecognized data\n')
    return
end
    
% Call appropriate routine
switch kind
    case 'data'
        switch modality
            case 'connectn'
                [NewFile] = store_data_connect(W, parent, comment);
            otherwise
                [NewFile] = store_data_other(W, parent, comment);
        end
    case 'result'
        switch modality
            case {'source','atlas'}
                [NewFile] = store_result_source(W, parent, comment);               
            case 'connectn'
                [NewFile] = store_result_connect(W, parent, comment, labels);
            case 'connecta'
                [NewFile] = store_result_connect(W, parent, comment, labels, iAtlas(1));
        end
end

% Disp
fprintf('\n\n****\tManual import : be_store_me ****\n****\tNew file in database:%s****\n****\tdone.\n\n',NewFile);
db_reload_studies(iStudy)

return



function [kind, modality, iS, iAtlas]   =   what_kind(W)

% STudy
[S,iS]  =   bst_get('Study');
Prot    =   bst_get('ProtocolInfo');
kind    =   [];
modality=   [];
iAtlas  =   [];

% is data?
if ~isempty(S.Channel)
    P   =   S.Channel.Modalities;
    Ch  =   load( be_fullfile(Prot.STUDIES, S.Channel.FileName), 'Channel' );
    if numel(Ch.Channel)==size(W,1)
        % Data matrix
        kind        =   1;
        modality    =   P{ii};
    end
    for ii  =   1 : numel(P)
        iCh =   find( strcmp({Ch.Channel.Type},P{ii}) );
        if numel(iCh)==size(W,1)
            % Data matrix
            kind        =   'data';
            modality    =   P{ii};
        elseif numel(iCh)==( sqrt(2*size(W,1)+1/4)+1/2 )
            % Data connectivity vector
            kind        =   'data';
            modality    =   'connectn';
        end
    end
end

% is result?
BSsub       =   bst_get('Subject', S.BrainStormSubject);
if ~isempty(BSsub)
    surfF   =   load( be_fullfile(Prot.SUBJECTS, BSsub.Surface(BSsub.iCortex).FileName ) );    
	nSC     =   cell2mat( cellfun(@(a) numel(a), {surfF.Atlas(:).Scouts}, 'uni', 0) );
    if size(surfF.Vertices,1)==size(W,1)
        % Results matrix
        kind        =   'result';
        modality    =   'source';
    elseif numel(surfF.Vertices,1)==( sqrt(2*size(W,1)+1/4)-1/2 )        
        % Results connectivity vector
        kind        =   'result';
        modality    =   'connectn';
    elseif any( nSC==size(W,1) )
        % Results matrix - atlas
        kind        =   'result';
        modality    =   'atlas';
        iAtlas      =   find( nSC==size(W,1) );
    elseif any( nSC==sqrt(2*size(W,1)+1/4)-1/2 )
        % Results connectivity - atlas
        kind        =   'result';
        modality    =   'connecta';
        iAtlas      =   find( nSC==( sqrt(2*size(W,1)+1/4)-1/2 ) );              
    end
end

return



function [NewFile] = store_result_source( W, parent, comment )

% Study
[S, iS] =   bst_get('Study');
iProt   =   bst_get('ProtocolInfo');
SUB     =   bst_get('Subject', S.BrainStormSubject);
AAT     =   load(be_fullfile(iProt.SUBJECTS,SUB.Surface(SUB.iCortex).FileName),'Atlas');

% Parent infos
if ~isempty( parent )
    INF =   load( be_fullfile(iProt.STUDIES, parent), 'Time' );   
else
   INF.Time = 1 : max(2,size(W,2));
end
if isempty(comment)
    comment     =   'Manual import (cf. be_store_me)';
end
if isempty(S.HeadModel)
    S.HeadModel(1).FileName = '';
end
if size(W,2)==1
    W = W*ones(1,2);
end

% Create a new signal from the extracted ridge lines
sMat.Time           =   INF.Time;
sMat.ImageGridAmp   =   W;
sMat.ImagingKernel  =   [];
sMat.Comment        =   comment;
sMat.nComponents    = 	1;
sMat.Function       =   'mem';
sMat.GridLoc        =   [];
sMat.Freqs          = 	'';
sMat.HeadModelFile  =   S.HeadModel.FileName;
sMat.HeadModelType  =   'surface';
sMat.DataFile       =   parent;
sMat                =   bst_history('add', sMat, 'compute', 'Data imported through be_store_me');
sMat.RefRowNames    =   1:size(sMat.ImageGridAmp,1);
sMat.RowNames       =   1:size(sMat.ImageGridAmp,1);
sMat.SurfaceFile    =   SUB.Surface(SUB.iCortex).FileName;

% File tag
fileTag = 'results_manual';

% ===== SAVE FILE =====

% Get output study
sOutputStudy    = bst_get('Study', iS);

% Output filename
NewFile         = bst_process('GetNewFilename', bst_fileparts(sOutputStudy.FileName), [fileTag]);
% Save file
save(NewFile, '-struct', 'sMat');
% Add file to database structure
db_add_data(iS, NewFile, sMat);

return



function [NewFile] = store_result_connect( W, parent, comment, labels, iAtlas )

% Study
[S, iS] =   bst_get('Study');
Prot    =   bst_get('ProtocolInfo');
SUB     =   bst_get('Subject', S.BrainStormSubject);
AAT     =   load(be_fullfile(Prot.SUBJECTS,SUB.Surface(SUB.iCortex).FileName),'Atlas');

if isempty(comment)
    comment = 'Connectn - computed externally';
end
if isempty(labels)
    labels = 1 : size(W,2);
end    

% ===== PREPARE OUTPUT STRUCTURE =====
% Create file structure
FileMat             =   db_template('timefreqmat');
FileMat.Comment     =   comment;
FileMat.DataType    =   'results';
FileMat.TimeBands   =   {'cohere' 0 1};
FileMat.Measure     =   'other';
FileMat.Method      =   'cohere';
FileMat.nAvg        =   1;
FileMat.TF          =   reshape( W, size(W,1), 1, size(W,2) );
FileMat.Freqs       =   labels;
FileMat.Time        =   [0 1];
FileMat.SurfaceFile =   SUB.Surface(SUB.iCortex).FileName;
FileMat             =   bst_history('add', FileMat, 'import', 'Connectivity map (manual import : be_store_me)');

% Atlas
if ~isempty(iAtlas)
    FileMat.Atlas       =   AAT.Atlas(iAtlas);
    FileMat.RefRowNames =   {AAT.Atlas(iAtlas).Scouts(:).Label};
    FileMat.RowNames    =   {AAT.Atlas(iAtlas).Scouts(:).Label};
end

% ===== OPTIMIZE STORAGE FOR SYMMETRIC MATRIX =====
FileMat.Options = struct;
% Determine if output matrix is symmetric
FileMat.Options.isSymmetric = 1;

% ===== SAVE FILE =====
% Output study, in case of average
FileMat.DataFile    =   parent;    

% File tag
fileTag = 'connectn_manual';

% Output filename
NewFile = bst_process('GetNewFilename', bst_fileparts(S.FileName), ['timefreq_' fileTag '_consensus']);
% Save file
save(NewFile, '-struct', 'FileMat');
% Add file to database structure
db_add_data(iS, NewFile, FileMat);

return