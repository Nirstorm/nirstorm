function [alpha, CLS, OPTIONS] = be_scores2alpha(SCR, CLS, OPTIONS, varargin)
%BE_SCORES2ALPHA computes the initial probability of a parcel being active in 
%   the MEM model using the MSP scores of the sources within each patch
%
%   INPUTS:
%       -   SCR     : vector of MSP scores with dimension Nsources
%       -   CLS     : vector of parcel labels for each source (1xNsources)
%       -   OPTIONS : 
%               model.alpha_method  :   method of intialization of the active           
%                   initial parcel active probabilities. (1=mean of MSP scores
%                   of sources within  parcel, 2=max, 3=median, 4=all 
%                   probabilities set to 0.5, 5=all prob. set to 1)
%
%               model.alpha_threshold:  threshold on the active probabilites. 
%                   All prob. < threshold are set to 0 (parcel not part of the 
%                   MEM solution
%
%   OUTPUTS:
%       -   ALPHA   : vector of probabilities (1xNparcels)
%       -   CLS     : cell array (1xNparcels). Each cel contains the indices of        
%                     the sources within that parcel
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

                                                        
alpha = zeros(size(SCR));
ALPHA_METHOD = OPTIONS.model.alpha_method;
alpha_threshold = OPTIONS.model.alpha_threshold;


% Progress bar
hmem = [];
if numel(varargin)
    hmem    = varargin{1}(1);
    st      = varargin{1}(2);
    dr      = varargin{1}(3);
end

for jj=1:size(SCR,2)
    clusters = CLS(:,jj);
    scores = SCR(:,jj);
    nb_clusters = max(clusters);
    curr_cls = 1;
    for ii = 1:nb_clusters
        idCLS = clusters==ii;
        switch ALPHA_METHOD
                
            case 1  % Method 1 (alpha = mean in the parcel)
                alpha(idCLS,jj) = mean((scores(idCLS)));
                
            case 2  % Method 2 (alpha = max in the parcel)
                 alpha(idCLS,jj) = max((scores(idCLS)));
                
            case 3  % Method 3 (alpha = median in the parcel)
                 alpha(idCLS,jj) = median((scores(idCLS)));
                
            case 4  % Method 4 (alpha = 1/2)
                 alpha(idCLS,jj) = 0.5;
                 
            case 5  % Method 4 (alpha = 1/2)
                 alpha(idCLS,jj) = 1;
                
            otherwise
                error('Wrong ALPHA Method')
        end
        
        % Check if cluster's active proba<threshold
        if alpha(idCLS,jj) < alpha_threshold
            CLS(idCLS,jj) = 0;
            alpha(idCLS,jj) = 0;
        else
            CLS(idCLS,jj) = curr_cls;
            curr_cls = curr_cls + 1;
        end
        
        % update progress bar
         if hmem
             prg = round( (st + dr * (jj - 1 + ii/nb_clusters) / (size(SCR,2))) * 100 );
             waitbar(prg/100, hmem, {[OPTIONS.automatic.InverseMethod, 'Step 1/2 : Running cortex parcellization ... ' num2str(prg) ' % done']});
         end
    end
    
end

% REMOVING CLUSTERS WITH ALPHA < APLHA_THRESHOLD
alpha(alpha > 0.8) = 1;


return



