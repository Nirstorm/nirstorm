function [clusters,time] = be_downsample_sources(OPTIONS,FILE)
% BE_DOWNSAMPLE_SOURCES: data-driven parcellisation of a sources intensity
% matrix. Preliminary step for computing connectivity measures
%
%% ==============================================
% Copyright (C) 2012 - LATIS Team
%
%  Authors: LATIS team, 2013
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

% Prepare variables
iP          =   bst_get('ProtocolInfo');
cortex      =   load( be_fullfile(iP.SUBJECTS,FILE.SurfaceFile), 'Faces', 'Vertices', 'VertConn');
nbS         =   size( cortex.Vertices,1 );
nbC         =   OPTIONS.sources.maxclus;
clusters    =   {};
    
% Step 1 : Find neighborhood 1
if ~isfield( cortex, 'VertConn' )
    cortex.VertConn     =   be_get_neighbor_matrix( nbS, cortex.Faces );
end
nanidx      =   ~any(isnan(real(FILE.ImageGridAmp)));
src_ampl    =   FILE.ImageGridAmp(:, nanidx);
neigh1_val  =   zeros(nbS,0,'single');
neigh1_idx  =   zeros(nbS,0,'single');
clear FILE
for ii      =   1 : nbS  
    n                           =   find( cortex.VertConn(1:nbS,ii) );
    neigh1_idx(ii,1:numel(n))   =   n;  
    switch OPTIONS.sources.metric
        case 'correlation'
            nei     =   corrcoef( src_ampl([ii; neigh1_idx(ii,1:numel(n))'],:)' );            
        case {'plv','iplv'}
            [d,nei] =   be_get_plv( src_ampl([ii; neigh1_idx(ii,1:numel(n))'],:), OPTIONS.sources.metric, 1 );            
    end
    neigh1_val( ii,1:numel(n) ) =   single( abs( nei(2:end,1) ) );     
end
clear cortex src_ampl

% Step 2 : Make linkage
Z   =   zeros( nbS-1, 3, 'single' );
nIx =   single( 1:nbS );
time=   [];
avL =   logical( ones(1,nbS) );
for ii = 1 : nbS-1
    
    % get max - light
    [tmp,tmpi]   =   max(neigh1_val(avL,:),[],2);
    maxima          =   zeros(nbS,1);
    imaxima         =   zeros(nbS,1);
    maxima(avL)     =   tmp;
    imaxima(avL)    =   tmpi;
    
    % if no link available, make one
    if ~sum(maxima)
        fprintf('New link\n')
        idnnul  =   find(nIx>0,2,'first');
        neigh1_idx( idnnul,: )    	=   0;
        neigh1_val( idnnul,: )     	=   0;
        neigh1_idx( idnnul(1),1 )   = idnnul(2);
        neigh1_idx( idnnul(2),1)    = idnnul(1);
        neigh1_val( idnnul(1),1 )   = 1/(10000*ii);
        neigh1_val( idnnul(2),1 )   = 1/(10000*ii);
        
        % get max - light
        avL(idnnul)	=   true;
        [tmp,tmpi] 	=   max(neigh1_val(avL,:),[],2);    
        maxima     	=   zeros(nbS,1);
        imaxima    	=   zeros(nbS,1);
        maxima(avL)	=   tmp;
        imaxima(avL)=   tmpi;
    end
    
    % find next pair to pool
    id_max1             =   find( maxima==max(maxima), 1, 'first' );
    id_max2_rel         =   imaxima( id_max1 );  
    id_max2             =   neigh1_idx(id_max1,id_max2_rel);
    id_max1_rel         =   imaxima(id_max2);
    if id_max1>id_max2
        tmp = [id_max1 id_max2_rel];
        id_max1     =   id_max2;
        id_max2     =   tmp(1);
        id_max2_rel =   id_max1_rel;
        id_max1_rel =   tmp(2);
    end
    
    % make cluster
    Z(ii,:)             =   [nIx(id_max1) nIx(id_max2) maxima(id_max1)];
    nIx(id_max1)        =   nbS + ii;
    nIx(id_max2)        =   0;
    neigh1_val(id_max1,id_max2_rel) = 0;
    
    %% update distances %%
    
    % - neighborhood
    nei1        =   neigh1_idx(id_max1,:);   
    nei2        =   neigh1_idx(id_max2,:);  
    neis        =   [nei1(~~nei1) nei2(~~nei2)];
    [dum,ix1]   =   unique(neis, 'first');
    neis        =   neis( sort(ix1) )';
    neis( neis==id_max1 ) = [];
    
    % - distance table
    dtable  =   zeros( 2, numel(neis) );
    in11    =   1 : sum(~~nei1);
    dtable(1, in11) =   neigh1_val( id_max1, in11 );
    [dum,in21,in22] =   intersect(nei2, neis);
    dtable(2, in22) =   neigh1_val( id_max2, in21 );
    dtable(3,:)     =   mean( dtable );
    
    % - unique cluster
    neigh1_idx( id_max1, 1:numel(neis) )    =   neis;
    neigh1_val( id_max1, 1:numel(neis) )    =   dtable(3,:);
    neigh1_val( id_max2, : )                =   0;
    avL(id_max2)                            =   false;
    
    % update neighbors distance
    if neis
        % find neighborhood with max1
        avaix   =   1:numel(neis);
        itable  =   neigh1_idx(neis,:);
        [r,c]   =   find( id_max1 == itable );
        if r
            [s,s]   =   sort(r);
            r       =   r(s);
            c       =   c(s);
            neigh1_val( (c-1)*nbS + neis(r) ) = dtable(3,r);
            avaix(r)=   [];
        end
        
        [r2,c2] =   find( id_max2 == itable(r,:) );
        if r2
            [s,s]   =   sort(r2);
            r2      =   r2(s);
            c2      =   c2(s);
            neigh1_val( (c2-1)*nbS + neis(r2) ) = 0;
        end
        
        % find neighborhood with only max2
        [r3,c3] =   find( id_max2 == itable(avaix,:) );
        if r3
            [s,s]   =   sort(r3);
            r3      =   avaix( r3(s) );
            c3      =   c3(s);
            neigh1_idx( (c3-1)*nbS + neis( r3 ) ) = id_max1;
            neigh1_val( (c3-1)*nbS + neis( r3 ) ) = dtable(3,r3);
        end
        
%         % make sure vectors are still sorted
%         minitable = neigh1_idx( neis, : )';
%         mindtable = neigh1_val( neis, : )';
%         minitable( minitable==0 ) = nbC + 1;
%         [minitabls, idx] = sort( minitable );
%         idx = idx + ones(size(idx,1),1)*[0:numel(neis)-1];
%         mindtable = mindtable(idx);
%         minitabls( minitabls==(nbC+1) ) = 0;
%         neigh1_idx( neis, : ) = minitabls';
%         neigh1_val( neis, : ) = mindtable';

    end    
 
end

% Make clustering
T = cluster( Z, 'maxclust', nbC );

% Put in appropriate format
for ii = 1 : max(T)
    clusters{ii} = find(T==ii);
end

return