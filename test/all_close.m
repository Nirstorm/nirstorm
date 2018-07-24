function flag = all_close(v1, v2, rtol, atol)

if nargin < 3
    rtol = 1e-5; %default relative tolerance
end

if nargin < 4
    atol = 1e-5; %default absolute tolerance
end

% Convert to vectors if needed:
if ndims(v1) == 2
    v1 = v1(:)';
end

if ndims(v2) == 2
    v2 = v2(:)';
end

flag = all(abs(v1 - v2) <= (atol + rtol * max(abs(v1), abs(v2))));

end