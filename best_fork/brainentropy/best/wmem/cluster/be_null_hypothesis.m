function active_mean = be_null_hypothesis(M,G,clusters,noise_std,varargin)
%BE_NULL_HYPOTHESIS initializes the mean of 'active state' law.
%   ACTIVE_MEAN = BE_NULL_HYPOTHESIS(M, G, CLUSTERS, NOISE_STD,MEANFLAG)
%   initializes the mean of the 'active state' law and returns it in 
%   ACTIVE_MEAN. The input parameter M is a column vector representing the 
%   mesures. G is the lead field matrix [NB_Sensors x Nb_dipoles].  
%   CLUSTERS is a 1xNB_CLUSTERS CELL. Each cell contains a [1xNB_dipoles] 
%   matrix that contains the index number of each dipole in the cluster.
%   NOISE_STD represents the standard deviation of the noise function.
%   MEANFLAG==true for initialization with the averaged null hypothesis, 
%   MEANFLAG==false for initialization with the null hypothesis
%
%   The method uses the methodoly described in :
%   Amblard, C., Lapalme, E., Lina, J. 2004, Biomagnetic Source Detection
%       by Maximyum Entropy and Graphical  Models, IEEE Transactions on 
%       Biomedical Engineering, vol. 51, no 3, p. 427-442.
%
%   The formulas are :
%       for a given clusker k:
%       if mean == true 
%       qk* = nGk' * ( sum_on_clusters(G * G') + 1/noise_std )^-1  * m
%       nGk = 1/nb_dipoles * sum(Gk)
%       note:
%           qk* is equal to the mean thus equal for each dipole
%       
%       if mean == false
%       qk* = Gk' * ( sum_on_clusters(G * G') + 1/noise_std )^-1  * m
%       note:  
%           qk* is not equal for each dipole
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


if ~isempty(varargin)
    meanflag=varargin{1};
else
    meanflag=true;
end
% Initialization
nb_clusters = max(clusters);
nb_sensors = size(M, 1);
active_mean = cell(1,nb_clusters);
inv_noise_std = diag(1./noise_std);


% SOLVING THE EQUATIONS

% Calculating the second term (the sum)
second_term = zeros(nb_sensors, nb_sensors);

for ii = 1:nb_clusters
    % extracting the lead field of the cluster
    G_cluster = G(:,clusters==ii);  

    % Calculating the sum
    second_term = second_term + (G_cluster * G_cluster');
end

% Inverting the second term (the sum) using 90% of the eigenvalues.
% Using this technique yields a more stable result than the builtin INV or
% PINV function of MATLAB.
[U,S,V] = svd(second_term + inv_noise_std);
eigen = diag(S);
I = find(cumsum(eigen.^2)./ sum(eigen.^2) > 0.9,1);
second_term = V(:,1:I) * diag(1./eigen(1:I)) * U(:,1:I)';

% Multiplying the second (the sum) term by the third term (M)
second_third_term = second_term * M;

% Completing the equation by multiplying the fist term with the last two.
% This multiplication will generate all the mu.
for ii = 1:nb_clusters    
    
    % extracting the lead field of the cluster
    G_cluster = G(:,clusters==ii);
    
    nb_dipoles_in_cluster = sum(clusters==ii);
    
    if meanflag == true
        % generating the first term
        nG_cluster = (1/nb_dipoles_in_cluster) * sum(G_cluster,2);    
    
        % Multiplication of the first term and the last two
        % The active mean is the same for all dipoles        
        g = ones(nb_dipoles_in_cluster,1) ...
            * nG_cluster' * second_third_term;          
    else
        % Multiplication of the first term and the last two
        % The active mean is not the same for all dipoles        
        g = G_cluster' * second_third_term;         
    end
   
    active_mean{ii} = g;
end


return

