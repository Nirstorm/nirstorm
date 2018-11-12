function varargout = process_nst_surfacedepth( varargin )

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
% Authors: ZhengChen Cai (2017-2018)



eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
% Description the process
sProcess.Comment     = 'Compute distances from cortex surface to head surface';
sProcess.Category    = 'Custom';
sProcess.SubGroup    = 'NIRS - wip';
sProcess.Index       = 1000;
sProcess.Description = '';
% Definition of the input accepted by this process
sProcess.InputTypes  = {'import'};
% Definition of the outputs of this process
sProcess.OutputTypes = {'import'};
sProcess.nInputs     = 1;
sProcess.nMinFiles   = 0;
sProcess.options.subjectname.Comment = 'Subject name:';
    sProcess.options.subjectname.Type    = 'subjectname';
    sProcess.options.subjectname.Value   = '';
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};
% Get subject name
if ~isempty(sProcess.options.subjectname.Value)
    SubjectName = file_standardize(sProcess.options.subjectname.Value);
else
    SubjectName = sInputs.SubjectName;
end
if isempty(SubjectName)
    bst_report('Error', sProcess, [], 'Subject name is empty.');
    return;
end

% ===== GET SUBJECT =====
% Get subject
[sSubject, iSubject] = bst_get('Subject', SubjectName);
if isempty(iSubject)
    bst_report('Error', sProcess, [], ['Subject "' SubjectName '" does not exist.']);
    return
end
% Save the new head model
    sStudy = sSubject;
    
    cortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    head = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);
    [~,dis]=knnsearch(head.Vertices,cortex.Vertices); 
    name = ['Distance from ' sSubject.Surface(sSubject.iCortex).Comment ' to ' sSubject.Surface(sSubject.iScalp).Comment ];
    iStudy = db_add_condition(sSubject.Name, 'distance');
    ResultFile = bst_process('GetNewFilename', bst_fileparts(sSubject.FileName), ...
                            ['/distance/results_' protect_fn_str(name)]);
                        
       
    ResultsMat = db_template('resultsmat');
    ResultsMat.Comment       = name;
    ResultsMat.Function      = '';
    ResultsMat.ImageGridAmp = dis;
    ResultsMat.Time          = 0;
    ResultsMat.DataFile      = [];%[bst_fileparts(sSubject.FileName),'/distance/results_' protect_fn_str(name) '.mat'];%//
    ResultsMat.HeadModelFile = []; 
    ResultsMat.HeadModelType = [];
    ResultsMat.ChannelFlag   = [];
    ResultsMat.GoodChannel   = [];
    ResultsMat.SurfaceFile   = sSubject.Surface(sSubject.iCortex).FileName;%//
    ResultsMat.GridLoc    = [];
    ResultsMat.GridOrient = [];
    ResultsMat.nAvg      = 1;
    % History
    ResultsMat = bst_history('add', ResultsMat, 'compute', 'distance from cortex to head surface calculated');
    % Save new file structure
    bst_save(ResultFile, ResultsMat, 'v6');
    sStudy_new = bst_get('Study', iStudy);
    
    % ===== REGISTER NEW FILE =====
    % Create new results structure
    newResult = db_template('results');
    newResult.Comment       = name;
    newResult.FileName      = file_short(ResultFile);
    newResult.DataFile      = ''; %sInputs.FileName;
    newResult.isLink        = 0;
    newResult.HeadModelType = 'surface';
    sStudy_new.Result = newResult;
    
    %iResult = length(sStudy.Result) + 1;
    
    %sStudy_new.Result = ResultsMat;
    % Update Brainstorm database
    bst_set('Study', iStudy, sStudy_new);
end

function sfn = protect_fn_str(s)
sfn = strrep(s, ' ', '_');
end