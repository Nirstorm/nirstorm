function [J,varargout] = be_jmne_NIRS(G,M,OPTIONS)
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

param1 = [0.1:0.1:1 1:5:100 100:100:1000]; 

display('cMEM, solving MNE by L-curve ... ');
p = OPTIONS.model.depth_weigth_MNE;
beta=0;
Sigma_s = diag(power(diag(G'*G)+beta,p)); % beta = 0 p= 0.5
W = sqrt(Sigma_s);
scale = trace(G*G')./trace(W'*W);       % Scale alpha using trace(G*G')./trace(W'*W)
alpha = param1.*scale;
Fit = [];
Prior = [];
for i = 1:length(param1)
    J = ((G'*G+alpha(i).*Sigma_s)^-1)*G'*M; % Weighted MNE solution
    Fit = [Fit,norm(M-G*J)];       % Define Fit as a function of alpha
    Prior = [Prior,norm(W*J)];          % Define Prior as a function of alpha
end
[~,Index] = min(Fit/max(Fit)+Prior/max(Prior));  % Find the optimal alpha
J = ((G'*G+alpha(Index).*Sigma_s)^-1)*G'*M;
if nargout > 1
    varargout{1} = param1;
end
disp('cMEM, solving MNE by L-curve ... done');
figure()
plot(Prior, Fit,'b.');
hold on;plot(Prior(Index), Fit(Index),'ro');
hold off
xlabel('Norm |WJ|');
ylabel('Residual |M-GJ|');
title('L-curve');
return
