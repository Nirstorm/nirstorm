function formats = nst_get_formats()
formats.pair_re = 'S(?<src_id>[0-9]{1,2})D(?<det_id>[0-9]{1,2})';
formats.pair_fmt = 'S%dD%d';
formats.chan_re = [formats.pair_re '(?<measure>WL\d+|HbO|HbR|HbT)'];
formats.pair_fmt = 'S%dD%d';
end