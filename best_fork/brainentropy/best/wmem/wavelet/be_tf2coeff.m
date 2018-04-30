function [coeffs,j,k] = time_freq2coeff( obj, time, frequency)
% TIME_FREQ2COEFF obtains coefficients from wavelet transform
%
%   INPUTS:
%       -   obj :   data strcuture
%       -   time:   time limits
%       -   frequency:   :   frequency range
%
%   OUTPUS:
%       -   coeffs: time-scale coefficients
%       -   j   :   scale
%       -   k   :   time
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


freq = obj.frequ_analyzed;
if frequency>freq(1)
    disp('---! requested frequency out of upper bound.');
    frequency = freq(1);
elseif frequency<freq(end)
    disp('---! requested frequency out of lower bound.');
    frequency = freq(end);
end
[I J] = sort(abs(frequency-freq));
j_found = J(1);

times = obj.time_analyzed;
for i=2:j_found
    times = (times(1:2:end)+times(2:2:end)) / 2;
end
if time>times(end)
    disp('---! requested time out of upper bound.');
    time = times(end);
elseif time<times(1)
    disp('---! requested time out of lower bound.');
    time = times(1);
end
[I J] = sort(abs(times-time));
k_found  = J(1);
coeffs_j = obj.coefficients('levels',j_found,'channel',[1:obj.nb_channels]);
coeffs   = coeffs_j(:,k_found);
j        = j_found;
k        = k_found;
end
