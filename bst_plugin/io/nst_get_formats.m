function formats = nst_get_formats()
formats.pair_re = 'S(?<src_id>\d+)D(?<det_id>\d+)';
formats.pair_fmt = 'S%dD%d';
formats.chan_re = [formats.pair_re '(?<measure>WL\d+|HbO|HbR|HbT)'];
formats.pair_fmt = 'S%dD%d';
end