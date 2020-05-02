function contrast = nst_make_basic_contrasts(events_order)

contrast = struct();
for ievent=1:length(events_order)
    contrast(ievent).label = events_order{ievent};
    con_vec = zeros(1, length(events_order));
    con_vec(ievent) = 1;
    contrast(ievent).vector = sprintf('[%s]', strip(sprintf('%d ', con_vec)));
end

end