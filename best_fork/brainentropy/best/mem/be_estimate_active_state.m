function [estimated_active_probability, estimated_amp] =  ...
    be_estimate_active_state(mem_struct,varargin)
%   [ESTIMATED_ACTIVE_PROBABILITY, ESTIMATED_ACTIVE_MEAN] = 
%   be_estimate_active_state(MEM_STRUCT) returns the estimation of the active 
%   probability of each cluster in ESTIMATED_ACTIVE_PROBABILITY and the
%   estimation of the mean of the activated state function in 
%   ESTIVATED_ACTIVE_MEAN. MEM_STRUCT is a structure with the
%   initial parameters of the MEM created with the function the function 
%   MEMSTRUCT.  Please see the help of this function for more 
%   information.  The field LAMBDA in MEM_STRUCT is the optimized lambda
%   which is considered LAMBDA*.
%
%   The formulas are :
%       note : * after a variable means estimated variable
%
%       q* = alpha* * mu + 
%           (1 - alpha*) * omega * G' * lambda +
%           alpha* * sigma * G' * lambda*
%
%       alpha* = alpha / (alpha + (1-alpha) * exp(F*)
%
%       F* = F0*(G' * lambda*) - F1(G' * lambda*)
%           where F0 is the inactive state and
%                 F1 is the active state
%       F0(xi) = 1/2 * xi' * omega * xi
%       F1(xi) = 1/2 * xi' * sigma * xi + xi' * mu
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


% Initialization
nb_clusters = numel(mem_struct.clusters);
estimated_active_probability = zeros(nb_clusters,1);
estimated_amp = cell(1, nb_clusters);

lambda = mem_struct.lambda;                    

% RESOLVING THE EQUATIONS FOR EVERY CLUSTERS
for ii = 1 : nb_clusters
    
    % variable change to reflect the notation in the reference paper
    alpha = mem_struct.clusters(ii).active_probability;
    mu = mem_struct.clusters(ii).active_mean;
    sigma = mem_struct.clusters(ii).active_var;
    omega = mem_struct.clusters(ii).inactive_var;
    G_cluster = mem_struct.clusters(ii).G;
    
    % Estimating F*
    xi = G_cluster' * lambda;
    xi_trans = xi';                              
    
    % F0 is set to a dirac by default (omega=0).
    if isempty(omega)
        F0=0;
    else
        F0 = 1/2 * xi_trans * omega * xi;             
    end
    F1 = 1/2 * xi_trans * sigma * xi + xi_trans * mu;
    
    F = F0 - F1;
    
    % Estimating alpha*
    estimated_alpha = alpha / (alpha + (1 - alpha) * exp(F));
    estimated_active_probability(ii) = estimated_alpha;
    
    % Estimating mu*
    %    estimated_mu = mu + ...
    %             ((1 - estimated_alpha) * omega * G_cluster' * lambda) ...
    %                  / estimated_alpha + ...
    %             sigma * G_cluster' * lambda
    %   estimated_active_mean{ii} = estimated_mu;
    
    
    % Estimating amplitudes*
    estimated_amp{ii} = estimated_alpha * mu + ...
        ((1 - estimated_alpha) * omega * xi) + ...
        estimated_alpha * sigma * xi;
end

return

