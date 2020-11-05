function s_pval = nst_format_pval(pval)

if pval < 0.001
    s_pval = sprintf('%1.1e', pval);
elseif pval < 0.01
    s_pval = sprintf('%1.3f', pval);
elseif pval < 0.1
    s_pval = sprintf('%1.2f', pval);
else
    s_pval = sprintf('%1.1f', pval);    
end

end