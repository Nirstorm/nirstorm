function [iDs, iDi] = get_ID( OPTIONS )
%GET_ID uses brainstorm built-in functions to get the study to which belongs the 
% current dataset
%
%   INPUTS:
%        -   OPTIONS
%
%   OUTPUTS:
%       -   iDs :   study iD
%       -   iDi :   study #
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



CND     = bst_fileparts(OPTIONS.optional.DataFile);
SUB     = bst_fileparts(CND);
% DEL     = numel(bst_fileparts(SUB)) + 2;
[dum, iDs] = bst_get( 'Study', be_fullfile(CND, 'brainstormstudy.mat') );

if nargout == 2
    global GlobalData
    sNm = [GlobalData.DataBase.ProtocolInfo(GlobalData.DataBase.iProtocol).STUDIES filesep];
    iDi = find( strcmpi( strrep(OPTIONS.optional.DataFile, sNm, ''), {dum.Data.FileName}) ); 
end
    
    
return
