function s = nst_misc_mirror_sig_bounds(s, nbp)
    s = [s((nbp:-1:1) + 1, :) ; s ; s((size(s,1)-1:-1:(size(s,1)-nbp)), :)]; 
end