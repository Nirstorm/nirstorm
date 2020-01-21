function [J, varargout] = solve_wmn(M, G , W, varargin)
%SOLVE_WMN solves the WMN inverse problem
%   J = SOLVE_WMN( M, G, W) solves the WMN inverse problem and
%   returns the activations in J. M is a column vector representing the 
%   mesures. G is the lead field matrix . W is the weighting parameter.  If
%   W = I this method will solve the minimum norm (MN) method.
%
%   More information on input/output parameters:
%       M (mesures):                    DOUBLE  [Nb Sensors x 1]
%       G (lead field):                 DOUBLE  [NB Sensors x Nb Dipoles]
%       W (wighting parameter:          DOUBLE  [Nb Dipoles x Nb Dipoles]
%
%       J (activation amplitudes):      DOUBLE  [Nb Dipoles x 1]
%
%   Formula used:
%
%       J = G' * (G * G')^-1 * M
%       J = W^-1*G' * (G * W^-1*G' + regularization)^-1 * M
%
%   Reference page in Help browser
%       <a href = "matlab:doc function_name">doc function_name</a>
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


%% INITIALIZATION
alpha = 0.01;
if numel(varargin)
    alpha = varargin{1};
end

%% SOLVING THE EQUATION

regularization = alpha .* eye(size(G,1));

inv_wGt = W^-1 * G'; 
pinvG = inv_wGt * (G * inv_wGt + regularization)^(-1);

J = pinvG * M;

%% END
return

%% HISTORY
%
% $ AUTHOR: LATIS
%   
% $ MARCH 5, 2009 Etienne Lemay
%   - First Version

