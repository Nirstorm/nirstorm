function [OPTIONS] = be_struct_copy_fields( OPTIONS, Def, fields, varargin )
% This functions fills missing option fields from DEF to OPTIONS without
% overriding existing fields.
% 
% INPUT :
% -------
%
%   OPTIONS     :   destination structure
%
%   Def         :   source structure
%
%
% OUTPUT :
% --------
%
%   OPTIONS     :   updated structure
%
%% ==============================================
% Copyright (C) 2011 - LATIS Team
%
% Authors: LATIS team, 2011
% revised: 03/2012
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


OPTIONS.InverseMethod   =   'MEM';

% set override
override        =   0;
if numel(varargin)>0
    override    =   varargin{1};
end

% accept call fields
if isempty( fields )
    fields      =   fieldnames(Def);
end

% Loop on all be-specific fields
for ii = 1 : numel(fields)
    if isfield(OPTIONS, fields{ii}) && isstruct( OPTIONS.(fields{ii}) )
        OPTIONS.(fields{ii}) =  do_copy_fields(OPTIONS.(fields{ii}), Def.(fields{ii}), override);
    else
        OPTIONS.(fields{ii}) =  Def.(fields{ii});
    end
end

return



function [sDest] = do_copy_fields(sDest, sSrc, override)


if (nargin < 3) || isempty(override)
    override = 1;
end

% No fields to add
if isempty(sSrc)
    return
% No fields in destination structure
elseif isempty(sDest)
    sDest = sSrc;
% Fields in both structures
else
    namesSrc = fieldnames(sSrc);
    for i = 1:length(namesSrc)
        if override || ~isfield(sDest, namesSrc{i}) || isempty(sDest.(namesSrc{i}))
            sDest.(namesSrc{i}) = sSrc.(namesSrc{i});
        end
    end
end


return