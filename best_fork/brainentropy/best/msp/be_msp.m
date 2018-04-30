function [scores, OPTIONS] = be_msp(M, Gstruct, OPTIONS)
%   BE_MSP returns a vector of MSP scores as proposed by Mattout(2005).
%   [OPTIONS, scores] = BE_MSP(M, Gstruct, OPTIONS) returns the MSP scores
%   of the sources associated, with the data M and the structure Gstruct. 
%   M is a matrix of dimension [number of sensors x number of time samples] 
%   Gstruct is a structure that contains the column-wise normalized G, 
%   the svd factors U and the diagoanal of singular values:
%   Gstruct.G is a matrix of dimension [number of sensors x number of 
%   sources]. 
%   OPTIONS contains a threshold used to filter the data and the forward 
%   operator (see be_memsolver_multiM and reference below).
%   
%   Output:
%   scores is a vector of dimension [#sources x 1].
%   
% 	Reference
%       Mattout, J., M. Pelegrini-Issac, L. Garnero et H. Benali. 2005. 
%       Multivariate source prelocalization (MSP): Use of functionally 
%       informed basis functions for better conditioning the MEG inverse 
%       problem. NeuroImage, vol. 26, no 2, p. 356-373.
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


% Default options settings
Def_OPTIONS.clustering = struct('MSP_R2_threshold', .92);

% Return empty OPTIONS structure
if (nargin == 0)
    OPTIONS = Def_OPTIONS;
    return
end
% Check field names of passed OPTIONS and fill missing ones with default values
OPTIONS = be_struct_copy_fields(OPTIONS, Def_OPTIONS, [], 0);
scores  =zeros(size(Gstruct.Gn,2),1);
clear Def_OPTIONS

% Check for NaN in the data
if sum( any(isnan(M)) )
   return
end

% Normalize the data matrix. Eliminate any resulting NaN.
Mn = bsxfun(@rdivide, M, sqrt(sum(M.^2, 1)));
Mn(isnan(Mn)) = 0;

% Project the normalized data on eigenvectors.
gamma = Gstruct.U'*Mn;

% Calculate the multiple correlation coefficients R2.
% NOTE:
%   This step is different than what is proposed in the reference.
%   R2 is defined as:
%   R2 = diag(gamma * pinv(gamma'*gamma) * gamma');
% Actually, gamma'*gamma is a large matrix close to the Identity. We thus
% neglect it. This is of no consequence on the final scores
R2 = diag(gamma*gamma');

% Reorder the singular values as a function of R2.
[temp, indices] = sort(R2,'descend');
lambda = Gstruct.lambda(indices);

% Select the columns of B as a function of the ordered singular values up
% to the threshold value.
i_T = indices( 1:find(cumsum(lambda)./sum(lambda ) >= OPTIONS.clustering.MSP_R2_threshold,1) );
Ut = Gstruct.U(:,sort(i_T));

% Create the projector.
Ms = Ut*Ut'*Mn;
Ps = Ms * pinv(Ms'*Ms) * Ms';
clear Ms Ut R2 gamma C Mn indices i_T M

% Calculate the MSP scores.
Ps2=Ps*Gstruct.Gn;
for i=1:size(Gstruct.Gn,2)
     scores(i)= Gstruct.Gn(:,i)'*Ps2(:,i);
end
return
