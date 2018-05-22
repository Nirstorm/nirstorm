function suffixes = get_installable_extras(manifest_dir, matlab_version)
% Append extra installation scenarios for specific versions of matlab.
% This prevents installing features relying on core matlab functions that
% are not available.
% It add extras corresponding to MANIFEST_only_from.<mat_release> where 
% <release_date> is same or older than given matlab.
%
% Args:
%    - manifest_dir (str):
%         Directory where to search for MANIFEST.*
%    [- matlab_version (str)]:
%         matlab version string up to which features can be installed.
%         Default is current matlab version
%
% Outputs:
%     - suffixes (cell of str):
%         List of extra scenario labels corresponding to
%         MANIFEST<_only_from...>  files that will install functions that 
%         are valid only from specific matlab versions.
%

if nargin < 2
    matlab_version = version();
end

items = dir(manifest_dir);
suffixes = {};
funcs_not_installed = {};
funcs_not_installed_requirements = {};
for iitem=1:length(items)
    item = items(iitem);
    [root, basename, ext] = fileparts(item.name);
    if ~item.isdir && strcmp(basename, 'MANIFEST_only_from')
        try
            manifest_mat_rdate = ext(2:end); % 2: to discard dot
            manifest_mat_version = rdate_to_version(manifest_mat_rdate);
            if  compare_matlab_versions(matlab_version, manifest_mat_version) >= 0
                suffixes = [suffixes ['_only_from.' manifest_mat_rdate]];
            else
                new_funcs = readlines(fullfile(manifest_dir, item.name));
                funcs_not_installed = [funcs_not_installed new_funcs];
                funcs_not_installed_requirements = [funcs_not_installed_requirements ...
                                                    repmat({manifest_mat_rdate}, ...
                                                           size(new_funcs,1), 1)];
            end
        catch ME
            if ~strcmp(ME.identifier, 'Nirstorm:BadMatlabVersionString')
                rethrow(ME);
            end
        end
    end
end
if ~isempty(funcs_not_installed)
    warning_msg = [];
    s_funcs_not_installed = '';
    for ifunc=1:length(funcs_not_installed)
        s_funcs_not_installed = [s_funcs_not_installed funcs_not_installed{ifunc} ...
                                 sprintf(' (requires >= %s)\n', funcs_not_installed_requirements{ifunc})];
    end
    warning('Current matlab version (%s) too old, these files won''t be installed:\n%s', ...
            matlab_version, s_funcs_not_installed);
    %TODO: warn about non installed functions!
end

end

function lines = readlines(fn)
f = fopen(fn);             
lines = textscan(f,'%s','delimiter','\n');
lines = lines{1};
fclose(f);
end