function [hp,hptab] = be_display_time_scale_boxes(obj,OPTIONS) 
%DISPLAY_ENTROPY_DROPS displays a graphical window with MEM entropy drops at 
% each  iteration
%
%   INPUTS:
%       -   obj
%       -   OPTIONS
%
%   OUTPUTS:
%      -    hp
%      -    hptab
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
% parametres: hp , Tmin, Tmax, J, N et boites
% hp : panel handle
% Tmin : time(0)
% Tmax : time(end)
% J : levels
% N : nomber of samples
% TFboxes is a structure .k ans .j are (j,k) coordinates



warning('off','all');
hp = figure;
hptab = uitabgroup; drawnow;
if ~strcmp(obj.data_type,'discrete_wavelet')
    disp('NO WAVELET DISPLAY');
    return
end
for ii = 1:length(OPTIONS.mandatory.DataTypes)
    onglet{ii} = uitab(hptab,'title',OPTIONS.mandatory.DataTypes{ii});
%    values = [];
%    values = sqrt(sum(obj.WData{ii}(:,OPTIONS.selected_samples(1,:)).^2,1)/size(obj.WData{ii},1));
%    values = (values-min(values))/(max(values)-min(values)); 

    Tmin = obj.t0 - (obj.info_extension.start-1)/OPTIONS.automatic.sampling_rate;
    Tmax = Tmin + (size(obj.data{ii},2)-1)/OPTIONS.automatic.sampling_rate;
    J    = size(OPTIONS.automatic.scales,2);
    N    = size(obj.data{ii},2);
    T  = Tmax-Tmin;
    e  = 0.05;

    
    hpc = uipanel('Parent', onglet{ii}, ...
                  'Units', 'Normalized', ...
                  'Position', [0.01 0.01 0.98 0.98], ...
                  'FontWeight','demi');

    set(hpc,'Title',[' Time-frequency boxes '],'FontSize',8);

    ax = axes('parent',hpc, ...
    'outerPosition',[0.01 0.01 0.98 0.98], ...
    'YTick',0.5:J-0.5, ...
    'YTickLabel',num2cell(1:J));

    axis(ax,'fill'); 
    box(ax,'on'); 
    xlim(ax,[Tmin Tmax]); 
    ylim(ax,[0 J]); 
    xlabel(ax,'time (s)'); 
    ylabel(ax,'scale j'); 
    hold on
    MMM = colormap(gray(size(OPTIONS.automatic.selected_values{ii},2)));

    selection=OPTIONS.automatic.Modality(ii).selected_jk;

    for b=1:length(OPTIONS.automatic.selected_values{ii})
    l = T/N*2^selection(2,b);
    rectangle('parent',ax, ...
              'Position',[Tmin+(selection(3,b)-1)*l,selection(2,b)-1+e,l,1-2*e],...
              'LineWidth',1, ...
              'LineStyle','-', ...
              'edgecolor','k', ...
              'FaceColor',MMM(b,:));
    pause(0.015)
    end
    hold off
    
    for scl = 1:length(OPTIONS.wavelet.selected_scales)
        titre = [' scale ' num2str(OPTIONS.wavelet.selected_scales(scl)) ' '];
        hhh = uitab(hptab,'title',titre);
        hpc = uipanel('Parent', hhh, ...
                  'Units', 'Normalized', ...
                  'Position', [0.01 0.01 0.98 0.98], ...
                  'FontWeight','demi');

        set(hpc,'Title',[' Wavelet coefficients (% of power) '],'FontSize',8);
        ax = axes('parent',hpc, 'outerPosition',[0.01 0.01 0.98 0.98]);
        sj = OPTIONS.wavelet.selected_scales(scl);
        bj = find(OPTIONS.automatic.Modality(ii).selected_jk(2,:)==sj);
        tt = OPTIONS.automatic.Modality(ii).selected_jk(6,bj);
        vv = OPTIONS.automatic.selected_values{ii}(2,bj);
        tv = OPTIONS.automatic.selected_values{ii}(3,bj);
        stem(tt,vv,'xk'); hold on; plot(tt(tv==1),vv(tv==1),'xr'); hold off
    end
    
    pause(3)
    
end
return