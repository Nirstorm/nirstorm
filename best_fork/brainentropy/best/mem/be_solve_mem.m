function [J, varargout] = be_solve_mem(mem_struct)
%SOLVE_MEM solves the MEM inverse problem
%   J = SOLVE_MEM(MEM_STRUCT) solves the MEM inverse problem and
%   returns the activations in J. MEM_STRUCT is a structure with the
%   initial parameters of the MEM created with the function the function 
%   MEMSTRUCT.  Please see the help of this function for more 
%   information.  
%
%   [J, RESULTS] = SOLVE_MEM(MEM_STRUCT) solves the MEM inverse problem and
%   returns the activations in J and a result structure in RESULTS.  The
%   structure RESULTS contains fields:
%       -iterations: number of iterations it took to minimize the function.
%       -lambda: the final lambda parameter after minimization.
%       -entropy: the entropy drop.
%       -active_probability: the final active probability of the clusters
%       (alpha)
%
%   The method uses the methodology described in :
%   Amblard, C., Lapalme, E., Lina, J. 2004. ' Biomagnetic Source Detection
%       by Maximum Entropy and Graphical  Models '. IEEE Transactions on 
%       Biomedical Engineering, vol. 51, no 3, p. 427-442.
%
% ==============================================
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


% no validation of input/output parameters for optimization.

% Initialization
nb_clusters = numel(mem_struct.clusters);
nb_sources = mem_struct.nb_sources;
J = zeros(1, nb_sources);

% MINIMIZING THE FREE ENERGY
% We save the OPTIMAL LAMBDA (LAMBDA*) in the MEM_STRUCT.

[mem_struct.lambda, entropy, iter] = be_minimize_free_energy(mem_struct);


% ESTIMATING ACTIVE PROBABILITY AND THE ACTIVE MEAN
% The active probability (alpha) and the active mean (mu) is estimated with
% the minimized LAMBDA obtained in the prevous section.
[active_probability, active_amp] = be_estimate_active_state(mem_struct);


% GENERATING AMPLITUDES
% multiplying the active mean times the active probability will generate
% the estimated amplitude
for ii = 1 : nb_clusters
    J(mem_struct.clusters(ii).indices) = active_amp{ii};
end

if nargout == 2
    results.iterations         = iter;
    results.lambda             = mem_struct.lambda;
    results.entropy            = entropy;
    results.active_probability = active_probability;
    varargout{1}               = results;
end

return


%