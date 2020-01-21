function OUT = be_bpfilter(IN, Sf, cut)
%WARNING : the band-pass filter cuts the frequencies from [0, cut(1)] U
% [cut(2) Sf/2]. Values in CUT are then filtered OUT.
%
%   INPUTS:
%       -   IN  : signal
%       -   Sf  : sampling frequency 
%       -   cut : cutoff frequency 
%
%   OUTPUTS:
%       -   OUT
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


    %% -------------------- CHECK INPUT ARGUMENTS ---------------------- %%
    OUT = IN;
    % Define cutoff frequencies
    if isempty(cut)
        return
    elseif ~isnumeric(cut)
        error('Cutoff frequency/ies must be numeric')
    elseif numel(cut) > 1
        cf1 = cut(1); cf2 = cut(2);
    elseif cut(1)==abs(cut(1))
        cf1 = 0; cf2 = cut(1);
    else
        cf1 = abs(cut(1)); cf2 = Sf/2;
    end
    
    % Checks cutoff frequencies
    if (cf1>Sf/2) || (cf2>Sf/2)
        error('Cutoff frequencies must be within nyquist range')
    end
    
    % Check spectral resolution
    if Sf/size(IN,2) >= diff([cf1 cf2])/2
        sprintf('In be_bpfilter : Insufficient time samples. Can''t filter data.')
        return
    end
    %%% ---------------  END OF CHECK INPUT ARGUMENTS ----------------- %%%
    
    
    
    %% ------------------------ FILTERING PART ------------------------- %%
    % Look for dsp toolbox
    Vs      =   ver;
    dsp     =   0;
    if any( strcmp('Signal Processing Toolbox', {Vs.Name}) )
        dsp =   1;
    end
    
    % with signal processing toolbox
    if dsp
        % Butterworth fiter design
        [b,a]   =   butter( 2, 2*[cf1 cf2]/Sf );
    
        % Apply Filter 
        OUT     =   filter(b,a,IN')';   
    
    else
        OUT     =   bst_bandpass(IN, Sf, cf1, cf2, 0, 1); 
        
    end
    %%% ------------------- END OF FILTERING PART --------------------- %%%
    

return