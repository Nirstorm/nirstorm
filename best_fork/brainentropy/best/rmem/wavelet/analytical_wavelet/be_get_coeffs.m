function [CFs] = be_get_coeffs(TF, frequencies, line, varargin)
%BE_GET_COEFFS extracts a ridge signal from a continuous time-frequency plane
%
%   INPUTS: 
%       -   TF      :   Time-frequency plane (cell, 1xNsensors)
%       -   frequencies:vector of frequencies of the TF plane
%       -   line    :   vector of coefficient indicies
%
%   OUTPUTS: 
%       -   CFs     :   complex coefficients
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

    % Set default
    [nbC, nbF, nbS] = size(TF);
    method      = 'strict';
    tolerance   = 2;
    if numel(varargin) && ~isempty(varargin{1})
        method = varargin{1};
    end
    if numel(varargin)==2 && ~isempty(varargin{2})
        tolerance = varargin{1};
    end
    
    switch method
        
        case 'tolerate'
            
            % Converts cell array into matrix form (if not already)
            if iscell(TF)
                TF = cells2matrix(TF);
            end
            
            % Checks input coherence
            if numel(size(TF)) ~= 3 || size(TF,2) ~= numel(frequencies)
                error(['Wrong inputs. Vector ''frequencies'' must match ' ...
                    'number of levels in ''TF'' variable.'])
            end
            
            % Prepare new signal matrix
            CFs = zeros( size(TF,1), numel(line) );
            
            % Determine min/max frequency tolerated
            curF = frequencies( mod(line, numel(frequencies)) );
            tosplit = find_closest( [curF-tolerance curF+tolerance], frequencies );
            maxF = tosplit( 1:numel(line) );
            minF = tosplit( numel(line)+1:end );
            init = ceil( line(1)/numel(frequencies) );
            
            % Fill the matrix
            for ii = 1 : numel(line)
                Frng = minF(ii) : maxF(ii);
                cand = TF(:, Frng, ii+init-1)';
                [dum, idx] = max( abs(cand), [], 1 );
                dec =  [0:size(CFs,1)-1] * numel(Frng);
                CFs(:,ii) = cand(idx+dec)';
            end
            

        case 'strict'           
            
            % Prepare new signal matrix
            CFs = zeros( 0, numel(line) );
            
            % Loop on modalities
            for hh = 1 : numel(TF)
                
                % Fill the matrix
                for ii = 1 : numel(TF{hh})
                    tf      =   TF{hh}{ii};
                    CFs     =   [CFs; tf(line)];
                end
            end
    end

return
