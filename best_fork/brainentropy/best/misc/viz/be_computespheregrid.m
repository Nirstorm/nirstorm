%% ===== SOURCE GRID =====
% THIS FUNCTION IS A SUBFONCTION of channel_extrapm.m
% IT WAS COPY/PASTED FOR USE WITH be_viz_topo.m
%
% Create source space for a given channel file
function sourcespace = ComputeSphereGrid(bfs_center, bfs_radius)
    % Build a full grid
    x = .01:.01:bfs_radius;
    x = [-x 0 x];
    [X,Y,Z] = meshgrid(x,x,x);
    % Keep only the points that are within the sphere
    D = sqrt((X.^2)+(Y.^2)+(Z.^2));
    i = find((D < bfs_radius) & (D > .03));
    sourcespace = [X(i) + bfs_center(1), Y(i) + bfs_center(2), Z(i) + bfs_center(3)];

% FRANCOIS
%     % Source locations
%     mne_src = tess_sphere(362);
%     mne_src = mne_src .* .5 * bfs_radius;
%     mne_src = bst_bsxfun(@plus, mne_src, bfs_center');
    
% SYLVAIN
%     % Number of sources to interpolate on
%     nsrc = length(srcChan);
%     [x,y,z] = sphere(ceil(sqrt(nsrc)));
%     % Center source locations and scale them to the best-fitting sphere found 
%     mne_src = [x(:),y(:),z(:)];
%     mne_src_orig = mne_src(unique(round(linspace(1,size(mne_src,1),nsrc))),:);
%     mne_src = mne_src_orig * .5 * bfs_radius;
%     nsrc = size(mne_src,1);
%     mne_src = mne_src + repmat(bfs_center',nsrc,1);
end