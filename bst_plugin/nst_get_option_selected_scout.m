function scout_selection = nst_get_option_selected_scout(options, label)
option_name = ['scout_sel_' label];
option_info = [option_name '_info'];
assert(isfield(options, option_name));
assert(isfield(options, option_info));
selection = options.(option_name).Value{1};
if ~isstruct(selection)
    if isempty(options.(option_info).Value)
        bst_error('Not scout selected');
    end
    scout_selection = options.(option_info).Value(selection);
else % a small hack to be able to directly define scout object 
     % when process is used in scripts
    scout_selection = selection;
end
end