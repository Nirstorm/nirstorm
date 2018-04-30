function [AVG] = be_avgcell(ARRAY, weights)
    
% Determine weights
if ~exist('weights', 'var') | isempty(weights)
    weights     =   ones(1, numel(ARRAY));
end
weights     =   weights / sum(weights);

AVG = sparse( size(ARRAY{1},1), size(ARRAY{1},2) );
nAR = numel(ARRAY);
for ii = 1 : nAR
    if size(ARRAY{1}) ~= size(ARRAY{ii})
        error('Cells contents have different dimensions')
    end
    AVG = AVG + ARRAY{ii} * weights(ii);    
end

return