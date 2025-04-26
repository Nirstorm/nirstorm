function s = nst_protect_fn_str(s)

    s = strrep(s, '|', '_');
    s = strrep(s, '"', '');
    s = strrep(s, ':', '_');
    s = strrep(s, '(', '_');
    s = strrep(s, ')', '_');
    s = strrep(s, '[', '_');
    s = strrep(s, ']', '_');
    s = strrep(s, '!', '_');
    s = strrep(s, '__','_');
    s = strrep(s, ' ', '_');
    
end
