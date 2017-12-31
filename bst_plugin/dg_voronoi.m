function [vor dist] = dg_voronoi(img, voxsize, seeds, distance, aniso)
% Geodesic Discrete Voronoi Diagram
% FORMAT [vor dist] = dg_voronoi(img, voxsize, seeds, distance, aniso)
%
% img       - binary image (2D or 3D) with :  > 0  : inside
%                                             <= 0 : outside 
% voxsize   - size of the voxels in each dimension 
%			  (default is [1.0 1.0 1.0])
% seeds     - n x 3 array with the positions of the n seeds
% distance  - type of distance to use
%             (default is 'd34')
% aniso     - anisotropic map
%
% vor   - Geodesic Discrete Voronoi diagram 
%         (label is equal to the position of the seed in 'seeds')
% dist  - Geodesic Distance map of img with seeds as objects
%_______________________________________________________________________
%
% Compute the geodesic discrete Voronoi Diagram of an image of labelled 
% objects using front propagation). The distance map is also available 
% on output.
%
% This is a cmex-compiled algorithm (dg_voronoi.c)
%_______________________________________________________________________
% @(#)dg_voronoi.m                Guillaume Flandin            01/10/11

error(sprintf('Missing MEX-file: %s', mfilename));
