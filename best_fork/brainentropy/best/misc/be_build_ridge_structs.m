function [sProcess, sInput] = build_ridge_structs(OPTIONS)
%BUILD_RIDGE_STRUCTS translates BEst ridge signal output to Brainstorm-compatible result file
%
%   INPUTS:
%        -   OPTIONS
%
%   OUTPUTS:
%       -   sProcess
%       -   sInput
%      
%% ==============================================   
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2011
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

% Make data structure
        sInput = struct;
        sInput.iStudy       =   get_ID(OPTIONS);
        [stdI, stdiD]       =   bst_get('Study', sInput.iStudy);
        iP                  =   bst_get('ProtocolInfo');
        fileN               =   strrep(OPTIONS.automatic.DataFile, iP.STUDIES, '');
        sInput.iItem        =   find(strcmp( fileN(2:end), {stdI.Data(:).FileName}) );
        sInput.FileName     =   stdI.Data(sInput.iItem).FileName;
        sInput.FileType     =   'data';
        sInput.Comment      =   stdI.Data(sInput.iItem).Comment;
        sInput.Condition    =   stdI.Condition{1};
        sInput.SubjectFile  =   stdI.BrainStormSubject;
        sInput.A            =   OPTIONS.optional.Data;
        sInput.TimeVector   =   OPTIONS.mandatory.DataTime;
        subI                =   bst_get('Subject', sInput.SubjectFile);
        chnI                =   bst_get('ChannelForStudy', stdiD);
        sInput.SubjectName  =   subI.Name;
        sInput.DataFile     =   '';
        sInput.ChannelFile  =   chnI.FileName;
        
        % Make process structure
        sProcess.options.timeW.Value        =   OPTIONS.optional.TimeSegment;
        sProcess.options.freqRNG.Value      =   {OPTIONS.ridges.frequency_range};
        sProcess.options.bsl.Value          =   OPTIONS.optional.BaselineSegment;
        sProcess.options.threshchx.Value    =   ~isempty(OPTIONS.optional.BaselineSegment);
        sProcess.options.fixed.Value        =   {OPTIONS.ridges.strength_threshold};
        sProcess.options.dur.Value          =   {OPTIONS.ridges.min_duration};
        
        n = numel(OPTIONS.mandatory.DataTypes);
        OPTIONS.mandatory.DataTypes(3:2:2*n)    =   OPTIONS.mandatory.DataTypes(2:end);
        OPTIONS.mandatory.DataTypes(2:2:2*n-1)  =   {' '};
        sProcess.options.sensortypes.Value      =   [OPTIONS.mandatory.DataTypes{:}];
        
return