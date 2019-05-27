function nst_bst_export_events(sFile, bst_events, event_fn)
% Wrapper around bst functions to load event in a backward compatible mode
% Load according to format based on extension

if ischar(bst_events)
    bst_evt_data = load(file_fullpath(bst_events), 'Events');
    bst_events = bst_evt_data.Events;
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