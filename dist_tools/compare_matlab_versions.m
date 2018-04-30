function order = compare_matlab_versions(v1, v2)
% Compare given matlab version strings.
% Return -1 if v1<v2, 0 if v1==v2, 1 if v1>v2
%
% Based on @see verLessThan

if ~ischar(v1) || ~ischar(v2)
    error(message('Nirstorm:compare_matlab_versions:invalidInput'))
end

v1_parts = getParts(v1);
v2_parts = getParts(v2);
for ipart=1:length(v1_parts)
    if v1_parts(ipart) ~= v2_parts(ipart)
        if v1_parts(ipart) < v2_parts(ipart)
            order = -1;
            return;
        else
            order = 1;
            return;
        end
    end
end
order = 0;
end

function parts = getParts(V)
parts = sscanf(V, '%d.%d.%d.%d')';
% zero-fills:
if length(parts) < 2
    parts(2) = 0;  
end
if length(parts) < 3
    parts(3) = 0; 
end
if length(parts) < 4
    parts(4) = 0;
end
end