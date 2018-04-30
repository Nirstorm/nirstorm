function [dwjk, OPTIONS] = be_wdenoise_csoft(wjk, info_extension, W, C, OPTIONS)
% soft schrinkage in the complex (or real) case for multivariate
% data (wjk) in the wavelet representation. We first whitten the data 
% (using W), then denoise, then restore the covariance (using C).
% Level per level, we consider the var-covar of the difference wjk-dwjk.
% This is stored in OPTIONS.wavelet_noise_variance{j}
%
%   INPUTS:
%       -   wjk :   discrete time-scale plane
%       -   W   :   whitening matrix
%       -   C   :   covariance matrix
%       -   OPTIONS:options structure (see be_main.m)
%
%   OUTPUTS:
%       -   dwjk:   denoised time-scale plane
%       -   OPTIONS
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: LATIS, 2012
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


% whittening and denoise:
[dwjk, OPTIONS] = be_denoise_csoft( W*wjk, info_extension, OPTIONS);
% restore the covariance:
dwjk = C*dwjk;
return