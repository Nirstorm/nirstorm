function ftable = nst_filter_table(table_in, filters)

table_cols = table_in.Properties.VariableNames;
if ismember('entry__', table_cols) % Are you that twisted??
    throw(MException('Nirstorm:UnsupportedTableWithEntry__', ...
                     '"Entry__" cannot be a table column'));
end

if isempty(fieldnames(filters))
    ftable = table_in;
    return;
end

check_filters(table_in, filters);

ftable = apply_filters(table_in, filters, 'include');
ftable = apply_filters(ftable, filters, 'exclude');

end

function ftable = apply_filters(table_in, filters, operation)
ftable = table_in(get_filters_selection(table_in, filters, operation),:);
end

function check_filters(data_table, filters)
table_cols = data_table.Properties.VariableNames;
allowed_filter_names = [table_cols 'entry__'];
filt_cnames = fieldnames(filters);
for ifilt=1:length(filt_cnames)
    cname = filt_cnames{ifilt};
    if ~any(strcmp(cname, allowed_filter_names)) 
        throw(MException('Nirstorm:InvalidFilterColumn', ...
                         ['Wrong filter "' cname ...
                          '". It must either be "entry" or a column name among: ' ...
                          strjoin(table_cols, ', ') '.']));
    end
    filter_def = filters.(cname);
    if ~isstruct(filter_def)
        throw(MException('Nirstorm:InvalidFilterType', ...
                         sprintf('Invalid filter definition for %s, must be struct with fields "include" or "exclude"', ...
                                 cname)));
    end
    filter_fields = fieldnames(filter_def);
    for ifield=1:length(filter_fields)
        operation = filter_fields{ifield};
        if ~any(strcmp(operation, {'include', 'exclude'}))  
               throw(MException('Nirstorm:InvalidFilterOperation', ...
                                ['Invalid filter "' cname '.' operation ...
                                 '". Allowed fields: "include" and "exclude"']));
        end
        % Ensure that filter is consistent with table specification
        filter = filter_def.(operation);
        if strcmp(cname, 'entry__')
            if ~isstruct(filter) || (isfield(filter, 'include') && ~all(ismember(fieldnames(filter.include), table_cols))) || ...
                (isfield(filter, 'exclude') && ~all(ismember(fieldnames(filter.exlude), table_cols)))
                throw(MException('Nirstorm:InvalidEntryFilter', ...
                                ['Invalid filter "' cname '.' operation ...
                                 '". Must be a struct with fields matching table columns']));
            end
        else
            if iscell(filter)
                if ~all(cellfun(@ischar, filter))
                    throw(MException('Nirstorm:InvalidFilterValue', ...
                                     ['Invalid filter "' cname '.' operation ...
                                      '". Must be a cell array of str']));
                end
                if ~iscell(data_table.(cname)) || ~ischar(data_table.(cname){1})
                    throw(MException('Nirstorm:InconsistentFilterType', ...
                                     ['Type of filter "' cname '.' operation ...
                                      '" not consistent with corresponding table column']));
                end
            elseif isnumeric(filter)
                 if ~isnumeric(data_table.(cname))
                    throw(MException('Nirstorm:InconsistentFilterType', ...
                                     ['Type of filter "' cname '.' operation ...
                                      '" not consistent with corresponding table column']));
                 end     
            elseif isa(filter, 'function_handle')
                 if isnumeric(data_table.(cname))
                     val1 = data_table.(cname)(1);
                 else
                     val1 = data_table.(cname){1};
                 end
                 filter_res1 = filter(val1);
                 if ~isscalar(filter_res1) || ~islogical(filter_res1)
                     throw(MException('Nirstorm:InvalidFilterFunction', ...
                         ['Invalid filter function for : ' cname '.' operation ...
                         '. It should be a predicate function (returning 0 or 1)']));
                 end
            else
                throw(MException('Nirstorm:InvalidFilter', ...
                    ['Bad filter for : ' cname '.' operation ...
                    '. It should be either a list of values to match or a predicate function']));
            end
        end
    end
end
end

function selection = get_filter_selection(data_table, cname, filter_def, operation)
selection = false(size(data_table, 1), 1);
filter = filter_def.(operation);
if strcmp(cname, 'entry__') %TOTEST
    entry_fields = fieldnames(filter);
    for ientry=1:length(filter)
        filter_entry = filter(ientry);
        subfilters = struct();
        for ifield=1:length(entry_fields)
            filter_entry_val = filter_entry.(entry_fields{ifield});
            if ischar(filter_entry_val)
                subfilters.(entry_fields{ifield}).include = {filter_entry_val};
            else
                subfilters.(entry_fields{ifield}).include = filter_entry_val;
            end
        end
        sel_entry = get_filters_selection(data_table, subfilters, 'include');
        nb_matches = sum(sel_entry);
        if nb_matches == 0 && strcmp(operation, 'include') %TOTEST
            throw(MException('Nirstorm:EntryNotFound', ...
                             ['Entry not found: ' newline evalc('disp(filter_entry)')]));
        elseif nb_matches > 1
            throw(MException('Nirstorm:NonUniqueEntry', ...
                             'Entry filter returned multiple selection'));
        end
        selection = selection | sel_entry;
    end
else
    if ~isa(filter, 'function_handle')
        if iscell(filter)
            for ival=1:length(filter)
                selection = selection | strcmp(data_table.(cname), filter{ival});
            end
        else
            for ival=1:length(filter)
                selection = selection | data_table.(cname) == filter(ival);
            end
        end
    else
        if iscell(data_table.(cname))
            selection = selection | cellfun(filter,data_table.(cname),'UniformOutput', 1);
        else
            selection = selection | arrayfun(filter,data_table.(cname));
        end
    end
end
end

function selection = get_filters_selection(acq_table, filters, operation)
filt_cnames = fieldnames(filters);

selection = true(size(acq_table, 1), 1);

%% Apply filters
for ifilt=1:length(filt_cnames)
    cname = filt_cnames{ifilt};
    filter_def = filters.(cname);
    if isfield(filter_def, operation)
        cur_selection = get_filter_selection(acq_table, cname, filter_def, ...
                                             operation);
        if strcmp(operation, 'include')
            selection = selection & cur_selection;
        else
            selection = selection & ~cur_selection;
        end
    end
end
end
