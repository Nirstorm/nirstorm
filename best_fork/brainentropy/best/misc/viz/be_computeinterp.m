%% ===== COMPUTE INTERPOLATION =====
% THIS FUNCTION IS A SUBFONCTION of channel_extrapm.m
% IT WAS COPY/PASTED FOR USE WITH be_viz_topo.m
%
% USAGE:  [WExtrap, src_xyz] = ComputeInterp(srcChan, destChan, bfs_center, bfs_radius, Whitener, epsilon)
%         [WExtrap, src_xyz] = ComputeInterp(srcChan, destChan, bfs_center, bfs_radius)
function [WExtrap, src_xyz] = be_computeinterp(srcChan, destChan, bfs_center, bfs_radius, Whitener, epsilon)
    % ===== PARSE INPUTS =====
    if (nargin < 6) || isempty(epsilon)
        epsilon = 0.0001;
    end
    if (nargin < 5) || isempty(Whitener)
        Whitener = [];
    end
    if (nargin < 4) || isempty(bfs_radius)
        bfs_radius = 0.07;
    end
    if (size(bfs_center,2) == 1)
        bfs_center = bfs_center';
    end
    WExtrap = [];

    % ===== COMPUTE SOURCE SPACE =====   
    % Computing spherical volume sourcespace.
    bfs_radius = 0.07;
    src_xyz = be_computespheregrid(bfs_center, bfs_radius)';
       
    % ===== COMPUTE LEADFIELDS =====
    % Single sphere model
    Param = struct('Center', bfs_center', 'Radii', bfs_radius);
    % Compute headmodels
    Gsrc2orig   = bst_meg_sph(src_xyz, srcChan,  repmat(Param, [1,length(srcChan)]));
    Gsrc2target = bst_meg_sph(src_xyz, destChan, repmat(Param, [1,length(destChan)]));
    if isempty(Gsrc2orig) || isempty(Gsrc2target)
        disp('EXTRAP> Error: One of the leadfields was not computed properly.');
        return
    end
    % Apply whitener to leadfield
    if ~isempty(Whitener)
        Gsrc2orig = Whitener * Gsrc2orig;
    end

    % Computing SVD of Grammian and truncating based on epsilon.
    [U,S,V] = svd(Gsrc2orig * Gsrc2orig');
    s = diag(S);
    ss = sum(s);
    sss = 1 - (cumsum(s)./ss);
    i = find(sss<epsilon);
    i = i(1);
    si = 1./s;
    si = si(1:i);
    %display(['Using ' num2str(i) ' singular values/vectors (i.e, dimensions)']);

    % Computing extrapolation matrix.
    coefs = V(:,1:i) * diag(si) * U(:,1:i)';
    WExtrap = Gsrc2target * Gsrc2orig' * coefs;

    % Apply whitener to the interpolator
    if ~isempty(Whitener)
        WExtrap = WExtrap * Whitener;
    end
end