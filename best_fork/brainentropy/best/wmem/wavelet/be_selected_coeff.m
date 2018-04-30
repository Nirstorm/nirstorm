function [OPTIONS] = be_selected_coeff(WM, obj, OPTIONS)
% BE_SELECTED_COEFF selects the wavelet coefficients to run the localization. 
% Only the boxes in the temporal interval of interest (TimeSegment) are kept
% (a box is kept if more than the half is included in the TimeSegment). 
% Only the boxes in the selected scale are kept.
% A maximum of 99% of the energy is kept.

% NOTES:
%     - This function is not optimized for stand-alone command calls.
%     - Please use the generic BST_SOURCEIMAGING function, or the GUI.a
%
% INPUTS:
%     - WM      : Data matrix (wavelet rep.) to be reconstructed using MEM
%     - obj     : Wavelet coefficients
%     - OPTIONS : Structure of parameters (described in be_main.m)
%
% OUTPUTS:
%     - Results : Structure
%          |- OPTIONS       : Keep track of parameters
%                           with the seleted wavelet coefficients
%
%
%% ==============================================   wavelet
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

% === Initialize the parameters
% number of scales (total)
nbl = size(OPTIONS.automatic.scales,2);
% number of time samples (extended)
nbb = size(WM{1},2);
% number of channels
nbc = size(WM{1},1);
% time t0 in the extended data
fs = OPTIONS.automatic.sampling_rate;
t0_extended = obj.t0 - (obj.info_extension.start -1)/fs;
% if OPTIONS.optional.verbose
%     disp(nbl)
%     disp(nbb)
%     disp(nbc)
%     disp(t0_extended)
% end

% === Table of selection
selected_jk = struct([]);
% we keep 99pc of power
pc_power = 0.99;

for jj = 1 : length(OPTIONS.mandatory.DataTypes)
    i_kept = [];
    j_kept = [];
    k_kept = [];
    t_kept = [];
    Wgfp   = [];
    selected_jk{jj} = [];
    for j=1:nbl
        deb    = 1+nbb/2^j + 2;
        fin    = nbb/2^(j-1)-2;
        Wgfp   = [Wgfp (sum(abs(WM{jj}(:,deb:fin)).^2,1)/nbc)];
        i_kept = [i_kept deb:fin];
        j_kept = [j_kept j*ones(1,length(deb:fin))];
        k_kept = [k_kept (deb-nbb/2^j):(fin-nbb/2^j)];
        t_kept = [t_kept t0_extended+(2^(j-1))*(2*((deb-nbb/2^j):(fin-nbb/2^j))-1)/fs ];
    end
    
    % === selection (based on the power)
        [Wgfp_sorted, I] = sort(Wgfp,'descend');
        Ic = find(cumsum(Wgfp_sorted)>=pc_power*sum(Wgfp_sorted),1,'first');
        i_kept = i_kept(I(1:Ic));
        j_kept = j_kept(I(1:Ic));
        k_kept = k_kept(I(1:Ic));
        t_kept = t_kept(I(1:Ic));
        e_kept = 100*Wgfp(I(1:Ic))/sum(Wgfp(I(1:Ic)));
        % what we finally keep:
        selected_jk{jj} = [i_kept ; j_kept ; k_kept ];
        
        % MSP windows (in the true data samples -not extended-)
        win_l = max((selected_jk{jj}(3,:)-1).*(2.^selected_jk{jj}(2,:))+1,obj.info_extension.start) - obj.info_extension.start +1;
        win_r = min( selected_jk{jj}(3,:).*(2.^selected_jk{jj}(2,:)),obj.info_extension.end)        - obj.info_extension.start +1;
    
    % === what we keep, finally:
    
        OPTIONS.automatic.Modality(jj).selected_jk = [selected_jk{jj} ; win_l ; win_r ; t_kept];
        selected_values{jj} = [Wgfp_sorted(1:Ic) ; e_kept];
end

% we keep the selection in the OPTIONS
OPTIONS.automatic.selected_values  = selected_values;

% Fusion of modalities in the selection: here we define OPTIONS.automatic.selected_samples
OPTIONS = be_fusion_of_samples(OPTIONS);

% we flag the boxes in the temporal interval of interest: 
t1 = OPTIONS.optional.TimeSegment(1  );
t2 = OPTIONS.optional.TimeSegment(end);
tmi = OPTIONS.automatic.selected_samples(6,:)-2.^(OPTIONS.automatic.selected_samples(2,:)-1)/fs/2;
tma = OPTIONS.automatic.selected_samples(6,:)+2.^(OPTIONS.automatic.selected_samples(2,:)-1)/fs/2;
sl = ~((tma<t1)|(tmi>t2));
OPTIONS.automatic.selected_samples = [OPTIONS.automatic.selected_samples ; sl];
for ii = 1 : length(OPTIONS.mandatory.DataTypes)
    OPTIONS.automatic.selected_values{ii} = [OPTIONS.automatic.selected_values{ii} ; sl];
end

% we flag the boxes with the scales of interest:
if ~isempty(OPTIONS.wavelet.selected_scales) && logical(prod(OPTIONS.wavelet.selected_scales))
    sl = ismember(OPTIONS.automatic.selected_samples(2,:),OPTIONS.wavelet.selected_scales);
else
    if OPTIONS.optional.verbose
         fprintf('%s, No specific scales selected\n', OPTIONS.mandatory.pipeline);
    end
    sl = ones(1,size(OPTIONS.automatic.selected_samples,2));
end
    OPTIONS.automatic.selected_samples = [OPTIONS.automatic.selected_samples ; sl];

% what do we keep finally ? coeff in the selected time window and selected
% scales:

sel = OPTIONS.automatic.selected_samples(8,:).*OPTIONS.automatic.selected_samples(9,:);
OPTIONS.automatic.selected_samples = OPTIONS.automatic.selected_samples(:,sel==1);
OPTIONS.automatic.selected_samples(8:9,:) = []; % (we can forget about the preselection, it is done)

% Keep only strongest t-f box
if OPTIONS.wavelet.single_box
    OPTIONS.automatic.selected_samples(:,2:end) = [];
end

if OPTIONS.optional.verbose
    fprintf('%s, wavelet selected boxes: you kept %d t-f boxes over a total of %d.\n', OPTIONS.mandatory.pipeline, size(OPTIONS.automatic.selected_samples,2), nbb); 
end
return