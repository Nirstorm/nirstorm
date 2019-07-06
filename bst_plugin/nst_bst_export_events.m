function nst_bst_export_events(sFile, event_fn)
% Wrapper around bst functions to export events.

if ischar(sFile)
    bst_evt_data = load(file_fullpath(sFile), 'Events');
    bst_events = bst_evt_data.Events;
else
    bst_events = sFile.Events;
end


[root, bfn, ext] = fileparts(event_fn);
switch ext
    case '.mat'
        s.events = bst_events;
        bst_save(event_fn, s, 'v7');
    otherwise
        error(sprintf('Unsupported extension for event file: "%s"', ext));
end
end