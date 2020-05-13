function nst_bst_save_events(bst_events, event_fn)
% Wrapper around bst functions to load event in a backward compatible mode
% Load according to format based on extension

[root, bfn, ext] = fileparts(event_fn);
switch ext
    case '.mat'
        wrapper.events = bst_events;
        % Emulate sampling frequency high enough to avoid rounding of times
        % wrapper.prop.sfreq = 2.0 / min(diff(sort([bst_events.times])));
        export_events(wrapper, [], event_fn);
    otherwise
        error(sprintf('Unsupported extension for event file: "%s"', ext));
end
end