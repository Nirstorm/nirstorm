function OPTIONS = be_fusion_of_samples(OPTIONS) 
% this function takes two modalities for which we computed the discete
% wavelet transform and selected coefficients of importance 
% (in OPTIONS.Modalitiy{}.selected_jk). It returns the selected samples for 
% the fusion of the two modalities. This will be stored in 
% OPTIONS.selected_samples that turns to be a table 3xnbr_of_sample. The 
% first raw is (as in the 1 modality case) the coeff index, the second and
% third raws are the position of of the selected coeff in mod{1} and mod{2}
% respectively. If the position is 0, the modality was absent.
% THE ORDER OF THE FINAL SELECTION IS RELATIVE TO THE MEAN POWER OF THE 
% OVERALL SELECTION
%   
%   INPUTS:
%       -   OPTIONS
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


if length(OPTIONS.mandatory.DataTypes)==2
    % in case of 2 modalities, the ordering will be related to the index of
    % the boxes (TO BE DONE)
    V1 = OPTIONS.automatic.Modality(1).selected_jk(4,:);
    V2 = OPTIONS.automatic.Modality(2).selected_jk(4,:);
    C1 = OPTIONS.automatic.Modality(1).selected_jk(1,:);
    C2 = OPTIONS.automatic.Modality(2).selected_jk(1,:);
    C12 = union(C1,C2);
    TAB = zeros(3,length(C12));
    VAL = zeros(2,length(C12));
    TAB(1,:) = C12;

    [tt,tt1,tt2]= intersect(C1,C12);
    TAB(2,tt2)  = tt1; 
    VAL(1,tt2) = V1(tt1);
    [tt,tt1,tt2]=intersect(C2,C12);
    TAB(3,tt2)  = tt1; 
    VAL(2,tt2) = V2(tt1);
    [v,i]=sort((VAL(1,:)+VAL(2,:))/2,'descend');
    OPTIONS.automatic.selected_samples = ...
        TAB(:,i);
else % only one modality, we keep the ordering related to the power 
    nboxes = size(OPTIONS.automatic.Modality(1).selected_jk,2);
    OPTIONS.automatic.selected_samples = ...
        [OPTIONS.automatic.Modality(1).selected_jk ; ones(1,nboxes)];
    % we may clear the OPTIONS.automatic.Modality(1).selected_jk
    % we code the modality in the last line of the table (1 here)
end
    
    
return