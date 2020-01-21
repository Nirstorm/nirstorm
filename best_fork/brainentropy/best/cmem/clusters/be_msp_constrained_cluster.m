function clusters = be_msp_constrained_cluster(Seeds,CPM,VertConn,OPTIONS)
%
%
% BE_MSP_CONSTRAINED_CLUSTER creates clusters associated to each seeds, 
% while preventing vertices associated to two different components
% to be associated with the same cluster (cf CPM matrix: defines the 
% constraints for what components could be fused and should not be fused)
%
% INPUTS: 
%        -Seeds                 : the vector containing the indices of the centers of gravity of each cluster
%        -CPM                   : Complementarity Probability Map is a vector of size equal to the number of vertices of the Mesh. 
%                                 This map ensures that vertices assigned
%                                 to distinct signal components are not assigned to the same parcel.
%        -VertConn              : Neighborhood matrix. 
%        -OPTIONS 				: Structure of parameters (described in be_cmem_solver.m)
% OUTPUTS: 
%        -clusters              : Variable whose each cell k contains the indices of the vertices belonging to each cluster k
%
%
%% ==============================================
% Copyright (C) 2011 - MultiFunkIm & LATIS Team
%
%  Authors: MultiFunkIm, LATIS team, 2011
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


% Neighbourhood matrix at a specific scale
%Adds connectivity on the diagonal (there should be zeros on the diagonal of nm
%initially)
VertConn=VertConn+speye(size(VertConn));


% Setting the neighbor matrix to the appropriate degree
neighborhood = VertConn;



if size(Seeds,1) > size(Seeds,2)
    Seeds = Seeds';
end


AlreadyAgregated = [];
AlreadyAgregated = unique([AlreadyAgregated,Seeds]);

clusters = cell(0,0);


%%%%%%%%%%% Initiation of the region growing loop  %%%%%%%%%%%%%%%%%
if OPTIONS.optional.verbose
	fprintf('%s, initialisation of the region growing loop ...',OPTIONS.mandatory.pipeline);
end

% 1- Seeds are already claissifed by decreasing order of priority

% 2- For each seed, on has to determine agregation variables 

for k = 1:length(Seeds)
    
    surrounding{k} = Seeds(k);
    clusters{k} = Seeds(k);       
        
end

if OPTIONS.optional.verbose
	fprintf(' done\n');
end

%%%%%%%%%%% Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% STOP CRITERIA of the LOOP
% - no more vertices to agregate 

% DEFINITIONS of the agregation criteria :
%
% - Homogeneity of a cluster  (thanks to the use of CPM map):
%       max( sum( CPM(ambiguous_dipole,cluster{k}) ) )
%
% - the SURROUNDING of  cluster k is composed of the closest neighnours minus the vertices 
% already agregated . As the \nconstruction of each SURROUNDING is done according to the 
% priority of each cluster, the frontier of the diffusion along the mesh is constrained 
% by criteria defined in the CPM matrix
%


it = 0;
STOP = 0;
if OPTIONS.optional.verbose
	fprintf('%s, we enter in the region growing loop\n',OPTIONS.mandatory.pipeline);
end
kk = 0;

while STOP == 0
        
    % 1- Update of the Agregation VariablesMise a jour des variables d'agregation
    it = it+1;              % Number of iteration of the region growing loop 
%     if Verbose == 1
%         disp(['be_msp_constrained_cluster:      - ITERATION No ',num2str(it)]);
%     end
    if OPTIONS.optional.verbose
	fprintf(' %d: ',it);
    end
    
    frontCom = [];          % Common frontiers  shared by all SURROUNDINGS at iteration it 
    ambiguous = [];         % ambiguous dipole that belongs to several SURROUNDINGS
    
    
    % 2- WE go over each cluster and we define the SURROUNDING and ambiguous dipoles    
%     if Verbose == 1
%         disp('be_msp_constrained_cluster:          --> definition of the surroundings and ambiguous dipoles...');
%     end
    if OPTIONS.optional.verbose
	fprintf('definition of the surroundings and ambiguous dipoles\n');
    end
    for k = 1:length(clusters)
        
        % 2.1- The surrounding of each cluster is composed of its
        % neighbours + previous SURROUNDING
        my_neighbours = [];
        for my_ind = 1:length(surrounding{k})
            
           my_neighbours = [my_neighbours unique( find(neighborhood( surrounding{k}(my_ind),:) ~=0) )];
        
        end
        my_neighbours = unique(my_neighbours);
        
        % 2.2- ...Vertices belonging to already formed clusters are being removed
        surrounding{k} = setdiff(my_neighbours,[0,clusters{k},AlreadyAgregated]);
        
        % 2.3- Definition of the common frontier and ambiguous points
        % NB: one cannot define non ambiguous points untill the common frontier is fully identified...
        
        isonfrontier = intersect(frontCom,surrounding{k});
        if isempty(isonfrontier) ~=1
            ambiguous = unique([ambiguous,isonfrontier]);
        end
        frontCom = unique([frontCom,surrounding{k}]);
        
    end
    
    
    
    %-----------------------------------------------%
    %- We agregate vertices belonging to the common frontier  -%    
    
    if isempty(frontCom) == 0           %---> it remains some vertices to agregate 
        
     % 3- Gestion of ambiguous points
        
        if isempty(ambiguous) ~= 1
            
%             if Verbose == 1
%                 disp('be_msp_constrained_cluster:           --> gestion of ambiguous vertices...'); 
%             end
            if OPTIONS.optional.verbose
                fprintf('    dealing with the ambiguous vertices\n');
            end
            
            
            class = zeros(length(ambiguous),length(clusters));
            tmp = 1:1:length(CPM);
            
            
            % 3.1- We �go over each cluster and we assess the CPM score of each ambiguous vertex
            
            for k = 1:length(clusters)
                
                
                list = surrounding{k};                                 % list of dipoles surrounding the cluster k 
                
                ambClus{k} = intersect(ambiguous,list);                 % liste of ambiguous vertice belonging to the surrounding of cluster k
                
                % Update the list for the current cluster 
                surrounding{k} = setdiff(list,ambClus{k});
                
                if isempty(surrounding{k}) == 0
                    
                    clusters{k} = [clusters{k},surrounding{k}];
                    AlreadyAgregated = [AlreadyAgregated,surrounding{k}];
                    
                end
                
                
                if isempty(ambClus{k}) == 0
                    
                    scoreK = sum(CPM(ambClus{k},clusters{k}),2)./length(clusters{k});          % score of cluster k
                    [c,ia,ib] = intersect(ambClus{k},ambiguous);                               % ib is estimated so that ambiguous(ib) = ambClusK;
                    class(ib,k) = scoreK;
                    
                end
                
            end
            
%             if Verbose == 1
%                 disp('be_msp_constrained_cluster:           -->   1. Esimation of the scores of the clusters: OK');
%             end
            if OPTIONS.optional.verbose
                fprintf('    estimation of the scores of the clusters\n');
            end
            
            
            % 4.2- Look for minima of the criterion to discriminate clusters 
            for k = 1:length(ambiguous)
                
                [tmp,i] = sort(-class(k,:));                                  % i contains the indices of the row k of the class ordered by increasing order
                stop2 = 0;
                l = 0;
                while stop2 == 0
                    
                    l = l+1;
                    if ismember(ambiguous(k),ambClus{i(l)}) == 1
                        
                        clusters{i(l)} = [clusters{i(l)},ambiguous(k)];
                        surrounding{i(l)} = [surrounding{i(l)},ambiguous(k)];
                        AlreadyAgregated = [AlreadyAgregated,ambiguous(k)];
                        
                        stop2 = 1;
                    end
                    
                end

            end
            
%             if Verbose == 1
%                 disp('be_msp_constrained_cluster:           -->  2. Looking for minima : OK');
%             end
            if OPTIONS.optional.verbose
                fprintf('    looking for minima (ok)\n');
            end
            
        else        % No more ambiguous vertices 
            

            if OPTIONS.optional.verbose
                fprintf('    trivial agregations\n');
            end
            
            for cl = 1:length(clusters)
                clusters{cl} = unique([clusters{cl},surrounding{cl}]);
                AlreadyAgregated = [AlreadyAgregated,surrounding{cl}];
            end
            
        end
        

        % 6- One check the conditions to exit the loop 
        if length(AlreadyAgregated) >= size(VertConn,1) 
            STOP = 1;
        end
        
        
    else                    %---> common frontier empty: no more points or seed pb ? 
        
        if length(AlreadyAgregated) == size(VertConn,1)
            STOP = 1;
        else    % Seed pb 
            if OPTIONS.optional.verbose
                fprintf('    agregation of remaining vertices\n');
            end
            %disp('be_msp_constrained_cluster:          --> agregations of the remaining vertices...')
            remaining = setdiff(1:size(VertConn,1),AlreadyAgregated);
            
            q = length(clusters);
            
            Seeds(q+1) = remaining(1);
            clusters{q+1} = remaining(1);
            surrounding{q+1} = remaining(1);  
            STOP2 = 0;
            
            while STOP2 == 0

                voi = setdiff(unique(neighborhood(clusters{q+1},:)),[0,AlreadyAgregated]);
                tmp = intersect(remaining,voi);
                if isempty(tmp) == 0
                    clusters{q+1} = [clusters{q+1},tmp];
                    AlreadyAgregated = [AlreadyAgregated,tmp];
                else
                    STOP2 = 1;
                end
                      
                
            end

            
            
        end
            
    end
    
   if it == 100
       if OPTIONS.optional.verbose
           fprintf('    too many iterations (WARNING)\n');
       end
       STOP = 1;
   end
    
end

if OPTIONS.optional.verbose
	fprintf('%s, Agregation of the clusters (ok)',OPTIONS.mandatory.pipeline);
end


