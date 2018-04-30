function [out] = be_fullfile(varargin)
% This function chains the strings in VARARGIN, separated by a file
% separator

if numel(varargin)>0
    out     =   varargin{1};
    for ii  =   2 : numel(varargin)
        out =   [out filesep varargin{ii}];
    end
end

return