function [param, OPTIONS] = be_msp_h0beta_fit(M, Gstruct, OPTIONS)
%   This function returns the parameters of the beta fit of the msp scores
%   computed from random windows in the baseline. The size of the windows
%   is set to 7 (3+1+3) and 20 windows are chosen.
%   Input:
%   G : lead fields matrix
%   M : data Nsensors x time from some baseline recording.
%   Output:
%   param.alpha and param.beta
%   
% 	Reference
%       JM Lina et al., IEEE 2011.
%   Authors: LATIS team, 2011.
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


% Default options settings
Def_OPTIONS = struct(...
    'Baseline_for_H0',          'true' ,...
    'FDR_max_window_size',      10, ...
    'FDR_min_window_size',      3, ...
    'FDR_window_number',        30);

% Return empty OPTIONS structure
if (nargin == 0)
    OPTIONS = Def_OPTIONS;
    return
end

% Check field names of passed OPTIONS and fill missing ones with default values
OPTIONS.clustering   =   be_struct_copy_fields(OPTIONS.clustering, Def_OPTIONS, [], 0);

%%

% Parameters:
NW    = OPTIONS.clustering.FDR_max_window_size;    % max size of the windows
nw    = OPTIONS.clustering.FDR_min_window_size;
NSAMP = OPTIONS.clustering.FDR_window_number;  % number of windows
% Initialisation:
n_win = unique(randi(size(M,2)-NW-1,1,3*NSAMP));
i_win = randperm(length(n_win));
i_win = i_win(1:min(NSAMP,length(n_win)));
n_win = n_win(i_win);
w_win = nw+floor((NW-nw)*rand(1,3*NSAMP));
w_win = w_win(i_win);
scores = zeros(size(Gstruct.Gn,2),1);

% averaged scores on the baseline's windows
for i = 1:length(i_win)
    scores = scores*(i-1);
    [scores_i, OPTIONS] = be_msp(M(:,n_win(i):n_win(i)+w_win(i)-1), Gstruct, OPTIONS);
    scores = (scores+scores_i)/i;
end
[param.alpha param.beta] = our_betafit(scores);
end

% ------ functions
function [alpha beta] = our_betafit(x)
% parameters of the beta fit
if ((min(x) <= 0) || (max(x) >= 1)) || (min(x) == max(x))
   disp('---! All msp values must be different, > 0 and < 1.');
   x = max(x,0);
   x = min(x,1);
end

n = length(x);
% Initial Estimates.
tmp1 = prod((1-x) .^ (1/n));
tmp2 = prod(x .^ (1/n));
tmp3 = (1 - tmp1 - tmp2);
ahat = 0.5*(1-tmp1) / tmp3;
bhat = 0.5*(1-tmp2) / tmp3;
pstart = [ahat bhat];

ld   = sum(log(x));
l1d  = sum(log(1-x));
phat = fminsearch(@(x) n*betaln(x(1),x(2))+(1-x(1))*ld+(1-x(2))*l1d, ...
       pstart,optimset('display','none'));
alpha = phat(1);
beta  = phat(2);
end
