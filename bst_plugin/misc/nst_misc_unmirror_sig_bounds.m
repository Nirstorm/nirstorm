function s = nst_misc_unmirror_sig_bounds(s, nbp)
    s = s((nbp+1):(size(s, 1)-nbp), :); 
end