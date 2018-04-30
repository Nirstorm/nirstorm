function OPTIONS = be_ridges_nullhypothesis(OPTIONS)
%BE_RIDGES_NULLHYPOTHESIS computes the distribution of the MSP scores of the 
% sources under the null hypothesis for a given ridge signal. This function 
% adapts parameters of BE_MODEL_OF_NULL_HYPOTHESIS.
%
%   INPUTS:
%        -   OPTIONS
%
%   OUTPUTS:
%       -   OPTIONS
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


if ~strcmpi(OPTIONS.mandatory.pipeline, 'rMEM') || ~isempty( OPTIONS.automatic.Modality(1).paramH0 )
    return
end

% Set variables
Frq     =   OPTIONS.wavelet.freqs_analyzed;
nbL     =   numel(Frq);
rdgFR   =   mean( Frq( mod( OPTIONS.automatic.selected_samples, nbL ) ) );
windS   =   OPTIONS.ridges.cycles_in_window / rdgFR * OPTIONS.automatic.sampling_rate;
OPTIONS.clustering.FDR_max_window_size = fix( 1.2 * windS ) + 1;
OPTIONS.clustering.FDR_min_window_size = fix( 0.8 * windS ) + 1;
OPTIONS =   be_model_of_null_hypothesis(OPTIONS);

return
