function [disl, groups] = be_localmaxima( slice, OPTIONS, varargin )

% VARS
[nbF,nbS]                   = size(slice);
lim1                        = 1;
lim2                        = nbF;
if numel(varargin)>0
    lim1    =   varargin{1}(1);
    lim2    =   varargin{1}(2);
end

% Retrive local maxima 
slicem1                     = [zeros(1,nbS); slice(1:end-1,:)];
slicep1                     = [slice(2:end,:); zeros(1,nbS)];
localmax                    = (slice>slicem1 & slice>slicep1);
clear slicem1 slicep1

% remove lowest freq
localmax(lim1,:)            = 0;
localmax(lim2,:)            = 0;

% Apply .95 threshold on local maxima at each time sample (see Amor,
% Rudrauf & al. 2005)
stdsl                     	= max( max(slice(lim1+1:lim2-1,:), [], 1) * OPTIONS.ridges.scalo_threshold , 1e-36);
stdsl                       = ones(nbF,1) * stdsl;
slice(~localmax)           	= 0;
disl                        = (slice>=stdsl);
clear stdsl slice

% Collect ridges lines - at least 3 samples long
minduration                 = fix(OPTIONS.ridges.min_duration*OPTIONS.automatic.sampling_rate/1000);
groups = [];
if any(disl(:))
    [groups]                = be_group_ridges(disl, minduration, OPTIONS.ridges.energy_threshold);
end

% Remove ridges lines w/less than 3 maxima
% idRemove                    = cell2mat( cellfun(@(a) numel(a)<3, groups, 'uni', 0) );
% disl( [groups{idRemove}] )  = [];
% disl                        = sparse(disl);
% groups(idRemove)            = [];

return