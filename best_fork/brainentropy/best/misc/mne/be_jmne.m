function [J,varargout] = be_jmne(G,M,OPTIONS)
% Compute the regularisation parameter based on what brainstorm already
% do. Note that this function replace be_solve_l_curve:
% BAYESEST2 solves the inverse problem by estimating the maximal posterior probability (MAP estimator).
%
%   INPUTS:
%       -   G       : matrice des lead-fields (donnee par le probleme direct)
%       -   M       : vecteur colonne contenant les donnees sur les capteurs
%       -   InvCovJ : inverse covariance matrix of the prior distribution
%       -   varargin{1} : param (alpha = param. trace(W*W')./trace(G*G')
%                   NB: sinon, alpha est evalue par la methode de la courbe en L
%
%   OUTPUTS:
%       -   J       : MAP estimator
%       -   varargout{1} : param
%       -   varargout{2} : pseudo-inverse of G
%% ==============================================
% Copyright (C) 2011 - Christophe Grova
%
%  Authors: Christophe Grova, 2011
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

[n_capt, n_sour] = size(G);

% selection of the data:
sample = be_closest(OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime);

M = M(:,sample(1):sample(2));

% We solve J = (W'W)^-1.G'.( G.(W'W)^-1.G' + alpha.Id )^-1.M
GG = G*G';
TrG   = trace(GG);

ratio = TrG/n_capt;

param1 = 1/9; %This parameter is the same as used in brainstorm
alpha = param1*ratio;

invG = G'*( GG + alpha.*eye(n_capt))^-1;
J = invG*M; %Regularisation parameter

if nargout > 1
    varargout{1} = param1;
end
return
