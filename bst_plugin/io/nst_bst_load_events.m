function  bst_events = nst_bst_load_events(event_fn)
% Wrapper around bst functions to save events in a backward compatible mode
% Save according to format based on extension
wrapper.events = [];
[root, bfn, ext] = fileparts(event_fn);
switch ext
    case '.mat'
        % Use dummy sampling frequency high enough so that times are not rounded
        wrapper.prop.sfreq = 1e6;
        wrapper = import_events(wrapper, [], event_fn, 'BST');
        bst_events = wrapper.events;
    otherwise
        error(sprintf('Unsupported extension for event file: "%s"', ext));
end
end