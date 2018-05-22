function joined = strjoin_(toks, delimiter)
% For compatibility
if ~verLessThan('matlab', '8.2')
    joined = strjoin(toks, delimiter);
else
    d = delimiter;
    n = numel(toks);
    if n == 0
        joined = '';
    elseif n == 1
        joined = toks{1};
    else
        ss = cell(1, 2 * n - 1);
        ss(1:2:end) = toks;
        [ss{2:2:end}] = deal(d);
        joined = [ss{:}];
    end
end
end