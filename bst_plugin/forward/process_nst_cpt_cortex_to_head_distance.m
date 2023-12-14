function varargout = process_nst_cpt_cortex_to_head_distance( varargin )
%  Compute the distance of each vertex of the cortical surface to the head
%  in mm

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
% Authors: Edouard Delaire (2023)

eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute distance from cortical surface to head surface';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS', 'Sources'};
    sProcess.Index       = 1406;
    sProcess.Description = '';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    % Definition of the outputs of this process
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 1;
    
    % Option: Subject name
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
if isfield(sProcess.options, 'subjectname') && ~isempty(sProcess.options.subjectname.Value)
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

% Obtain the cortical surface
sCortex = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
sHead   = in_tess_bst(sSubject.Surface(sSubject.iScalp).FileName);

distance = Compute(sCortex, sHead);

iStudy = db_add_condition(SubjectName, 'Distance');
sStudy = bst_get('Study', iStudy);


OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), ...
                        ['results_distance_cortex_to_head']);
                    
% ===== CREATE FILE STRUCTURE =====
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = 'Cortex to Head distance';
ResultsMat.Function      = '';
ResultsMat.Time          = [0];
ResultsMat.ImageGridAmp = distance*1000;
ResultsMat.ChannelFlag   = [];
ResultsMat.GoodChannel   = [];
ResultsMat.DisplayUnits  = 'mm';
ResultsMat.SurfaceFile   = sSubject.Surface(sSubject.iCortex).FileName;
% History
ResultsMat = bst_history('add', ResultsMat, 'compute', 'Compute distance fron cortex to head');
% Save new file structure
bst_save(OutputFile, ResultsMat, 'v6');
% Update database
db_add_data(iStudy, OutputFile, ResultsMat);
OutputFiles = {OutputFile};

end

function distance = Compute(surfaceA, surfaceB)
% For every point i in surfaceA, distance(i) is the minimum distance from
% that point to anypoint in surfaceB

    x = surfaceA.Vertices;
    y = surfaceB.Vertices;
    d = nst_pdist(x,y);
    distance = min(d,[],2);
end

