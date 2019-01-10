function item = nst_parse_bst_item_name(item_name)

%TODO: handle error scenarios

item.comment = '';
item.subject_name = '';
item.condition = '';

split_re = '(?<root>.*/)?(?<item>[^/]+)';
toks = regexp(item_name, split_re, 'names');

if ~isempty(toks)
    item.comment = toks.item;
    if ~isempty(toks.root)
        toks = regexp(toks.root(1:end-1), split_re, 'names');
        item.subject_name = toks.root(1:end-1);
        item.condition = toks.item;
    end
end

end