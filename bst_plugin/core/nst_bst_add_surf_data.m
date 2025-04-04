function [sStudy, ResultFile] = nst_bst_add_surf_data(data, time, head_model, file_tag, comment, ...
                                                      sInputs, sStudy, history_comment, surface_file, ...
                                                      sparse_storage, extra)
                                                  
                                                  
if nargin < 9 || isempty(surface_file)
    if ~isempty(head_model)
        surface_file =  file_short(head_model.SurfaceFile);
    else
        error('Surface file not defined');
    end
end

if nargin < 10
    sparse_storage = 0;
end

if nargin < 11 
    extra = struct();
end

%TODO: check consistency between data and nb of vertices

%% Save a cortical map to brainstorm with given data

ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                         ['results_' protect_fn_str(file_tag)]);

% ===== CREATE FILE STRUCTURE =====
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = comment;
ResultsMat.Function      = '';
ResultsMat.ImageGridAmp  = data;
ResultsMat.Time          = time;


if ~isempty(sInputs) 
    if isfield(sInputs, 'DataFile')
        ResultsMat.DataFile = sInputs.DataFile;
    elseif isfield(sInputs, 'FileName')
        ResultsMat.DataFile  = sInputs.FileName;
    end
end

if ~isempty(head_model)
    if ~isempty(sStudy.iHeadModel)
        ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    end
    ResultsMat.HeadModelType = head_model.HeadModelType;
end


ResultsMat.ChannelFlag   = [];
ResultsMat.GoodChannel   = [];
ResultsMat.SurfaceFile   = surface_file;
ResultsMat.GridLoc       = [];
ResultsMat.GridOrient    = [];
ResultsMat.nAvg          = 1;

% Add extra fields
 extra_fields = fieldnames(extra);
 for ifield = 1:length(extra_fields)
%      assert(~isfield(ResultsMat, extra_fields{ifield}));
     ResultsMat.(extra_fields{ifield}) = extra.(extra_fields{ifield});
 end

% History
if ~isempty(sInputs)
    sData = in_bst_data(sInputs.FileName,'History');
    if isfield(sData,'History')
        ResultsMat.History = sData.History;
    end
end

ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
% Save new file structure
bst_save(ResultFile, ResultsMat, 'v6');
% ===== REGISTER NEW FILE =====
% Create new results structure
newResult = db_template('results');
newResult.Comment       = comment;
newResult.FileName      = file_short(ResultFile);
newResult.isLink        = 0;
newResult.HeadModelType = ResultsMat.HeadModelType;
% Add new entry to the database
iResult = length(sStudy.Result) + 1;
sStudy.Result(iResult) = newResult;
% Update Brainstorm database
if ~isempty(sInputs) && isfield(sInputs,'iStudy')
    bst_set('Study', sInputs.iStudy, sStudy);
else
    [tmp, iStudy] = bst_get('Study', sStudy.FileName);
    bst_set('Study', iStudy, sStudy);
end
                                                  
end

function sfn = protect_fn_str(s)
sfn = strrep(s, ' | ', '--');
sfn = strrep(s, ' : ', '--');
sfn = strrep(s, ' :', '--');
sfn = strrep(s, ' ', '_');
end