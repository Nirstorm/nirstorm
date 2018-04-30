function [OPTIONS] = be_remove_dc(OPTIONS)
% This function removes the DC offset from the data using the baseline
% segment defined in the panel
%
%% ==============================================
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2014
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


% Compute mean during baseline
mu      =   mean( OPTIONS.optional.Baseline,2 );

% subtract from baseline
nSb     =   size( OPTIONS.optional.Baseline,2 );
muM     =   mu * ones(1, nSb);
OPTIONS.optional.Baseline = OPTIONS.optional.Baseline - muM;

% subtract from data
nSd     =   size( OPTIONS.mandatory.Data,2 );
muM     =   mu * ones(1, nSd);
OPTIONS.mandatory.Data = OPTIONS.mandatory.Data - muM;

return