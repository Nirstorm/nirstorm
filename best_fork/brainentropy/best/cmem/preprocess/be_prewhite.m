function [OPTIONS]  =  be_prewhite(OPTIONS)
% BE_PREWHITE prewhitens the data using empty room recordings
%
%   INPUTS:
%       -OPTIONS     : Structure of parameters (described in be_main.m)
%
%   OUTPUTS:
%       - OPTIONS	 : Updated structure
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

% 3 options : nothing, from the baseline, from the empty-room data
% if nothing: we have to specify the cov mat in the MEM solver
% if ~nothing: the covariance is identity in MEM

for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
   
    automatic.Modality(ii).units.Gain_units =  1;   
    automatic.Modality(ii).units.Data_units =  1;
    automatic.Modality(ii).units.Cov_units =  1; 
end

% Try empty room normalization
done = 0;
for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
    automatic       =   OPTIONS.automatic;
    [isdone, ERdata]=   check_emptyroom(OPTIONS, ii);
    done            =   min(done,isdone);
    whiteMat        =   diag( diag(cov(ERdata')).^(-1/2) );
    automatic.Modality(ii).gain             =   whiteMat * automatic.Modality(ii).gain;
    automatic.Modality(ii).units.Gain_units =   whiteMat;   
    automatic.Modality(ii).units.Data_units =   whiteMat;
    automatic.Modality(ii).units.Data_units =   whiteMat.^2; 
end

return    


function [done,ERdata]  =   check_emptyroom(OPTIONS, ii)

done = 1;
if size( OPTIONS.optional.EmptyRoom_data,1 )== numel(OPTIONS.mandatory.ChannelTypes)
    % No need for EmptyRoom_channels
    iChan   =   find( strcmp(OPTIONS.mandatory.DataTypes{ii}, OPTIONS.mandatory.ChannelTypes) );
    ERdata  =   OPTIONS.optional.EmptyRoom_data(iChan,:);
elseif ~isempty(OPTIONS.optional.Channel)
    % Check in EmptyRoom_channels for a match with data    
    chN     =   {OPTIONS.optional.Channel.Name};
    [a,b,c] =   intersect(OPTIONS.optional.EmptyRoom_channels, chN);
    [d,e]   =   sort(b);
    ERdata  =   OPTIONS.optional.EmptyRoom_data( c(e),: );
else
    % No usable empty room data
    ERdata  =   1;
    done    =   0;
end

return