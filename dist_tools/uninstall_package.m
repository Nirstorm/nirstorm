function uninstall_package(package_name, target_dir)
uninstall_script_fns = find_uninstall_scripts(package_name, target_dir);
if ~isempty(uninstall_script_fns)
    cellfun(@(uscript) run(uscript), uninstall_script_fns);
    cellfun(@(uscript) delete(uscript), uninstall_script_fns);
end
end

function found_script_fns = find_uninstall_scripts(package_name, folder)
found_items = dir(fullfile(folder, ['uninstall_' package_name '*.m']));
found_script_fns = {};
for iitem=1:length(found_items)
    if ~found_items(iitem).isdir
        found_script_fns{end+1} = fullfile(folder, found_items(iitem).name);
    end
    clear(found_items(iitem).name);
end
end