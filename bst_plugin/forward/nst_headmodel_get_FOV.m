function valid_nodes = nst_headmodel_get_FOV(ChannelMat, cortex, thresh_dis2cortex, ChannelFlag)
    %% define the reconstruction FOV
    
    if nargin < 4
        ChannelFlag = ones(1, length(ChannelMat.Channel));
    end
    
    montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel, ChannelFlag);
    src_locs = montage_info.src_pos;
    det_locs = montage_info.det_pos;
    
    optodes_pos = [src_locs;det_locs];
    % inflate surface 100% to calculate distances to optodes (see BST folder figure_3d.m line 2595)
    iVertices = 1:length(cortex.Vertices);
    % Smoothing factor
    SurfSmoothIterations = ceil(300 * 1 * length(iVertices) / 100000);
    % Calculate smoothed vertices locations
    Vertices_sm = cortex.Vertices;
    Vertices_sm(iVertices,:) = tess_smooth(cortex.Vertices(iVertices,:), 1, SurfSmoothIterations, cortex.VertConn(iVertices,iVertices), 1);
    dis2cortex = pdist2(Vertices_sm,optodes_pos);

    valid_nodes = find(min(dis2cortex,[],2)<thresh_dis2cortex);
end