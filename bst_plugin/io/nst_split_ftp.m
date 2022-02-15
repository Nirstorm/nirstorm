function [ftp_site, rfn] = nst_split_ftp(ftp)

ftp = strrep(ftp, 'ftp://', '');
toks = strsplit(ftp, '/');
ftp_site = toks{1};
if length(toks) > 1
    rfn = strjoin(toks(2:end), '/');
else
    rfn = '';
end

end