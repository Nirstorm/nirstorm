function data = webread(url)
% Alternative version of webread for matlab < 2014b
% Fall back to urlread
data = urlread(url);
end