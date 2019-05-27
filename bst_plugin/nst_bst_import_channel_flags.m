function nst_bst_import_channel_flags(sFile,  channel_flags_fn)
%NST_BST_EXPORT_CHANNEL_FLAGS Save channel flags from given sFile 
%                             into given channel_flags_fn.
% TODO: also save list of channels
if ischar(sFile)
    bst_chan_data = load(file_fullpath(sFile), 'ChannelFlag');
    bst_chan_flags = bst_chan_data.ChannelFlag;
end
[root, bfn, ext] = fileparts(channel_flags_fn);
switch ext
    case '.mat'
        s.channel_flags = bst_chan_flags;
        bst_save(channel_flags_fn, s, 'v7');
    otherwise
        error(sprintf('Unsupported extension for channel file: "%s"', ext));
end
end


