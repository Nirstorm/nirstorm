function [sStudy, ResultFile] = nst_bst_add_surf_data(data, time, head_model, name, ...
                                                      sInputs, sStudy, history_comment, surface_file)
                                                  
                                                  
if nargin < 8
    if ~isempty(head_model)
        surface_file =  file_short(head_model.SurfaceFile);
    else
        error('Surface file not defined');
    end
end

%% Save a cortical map to brainstorm with given data

ResultFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
    ['results_' protect_fn_str(name)]);

% ===== CREATE FILE STRUCTURE =====
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = name;
ResultsMat.Function      = '';
ResultsMat.ImageGridAmp = data;
ResultsMat.Time          = time;
ResultsMat.DataFile      = sInputs.FileName;
if ~isempty(head_model)
    if ~isempty(sStudy.iHeadModel)
        ResultsMat.HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    end
    ResultsMat.HeadModelType = head_model.HeadModelType;
else
    
end
ResultsMat.ChannelFlag   = [];
ResultsMat.GoodChannel   = [];
ResultsMat.SurfaceFile   = surface_file;
ResultsMat.GridLoc    = [];
ResultsMat.GridOrient = [];
ResultsMat.nAvg      = 1;
% History
ResultsMat = bst_history('add', ResultsMat, 'compute', history_comment);
% Save new file structure
bst_save(ResultFile, ResultsMat, 'v6');
% ===== REGISTER NEW FILE =====
% Create new results structure
newResult = db_template('results');
newResult.Comment       = name;
newResult.FileName      = file_short(ResultFile);
newResult.DataFile      = ''; %sInputs.FileName;
newResult.isLink        = 0;
newResult.HeadModelType = ResultsMat.HeadModelType;
% Add new entry to the database
iResult = length(sStudy.Result) + 1;
sStudy.Result(iResult) = newResult;
% Update Brainstorm database
bst_set('Study', sInputs.iStudy, sStudy);
                                                  
end

function sfn = protect_fn_str(s)
sfn = strrep(s, ' ', '_');
end