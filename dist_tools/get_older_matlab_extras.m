function suffixes = get_older_matlab_extras(manifest_dir, matlab_version)
% Append extra installation scenarios for old versions of matlab. These
% scenarios should provide alternate implementation of matlab core functions 
% not available for the given version.
% Add extras corresponding to MANIFEST_compat.<mat_release> where <release_date>
% is more recent than the given matlab version.
%
% Args:
%    - manifest_dir (str):
%         Directory where to search for MANIFEST.*
%    [- matlab_version (str)]:
%         matlab version string that must be supported.
%         Default is current matlab version
%
% Outputs:
%     - suffixes (cell of str):
%         List of extra scenario labels corresponding to MANIFEST<_suffix>
%         files that will install functions to support an older matlab
%

if nargin < 2
    matlab_version = version();
end

items = dir(manifest_dir);
suffixes = {};
for iitem=1:length(items)
    item = items(iitem);
    [root, basename, ext] = fileparts(item.name);
    if ~item.isdir && strcmp(basename, 'MANIFEST_compat')
        try
            manifest_mat_rdate = ext(2:end); % 2: to discard dot
            manifest_mat_version = rdate_to_version(manifest_mat_rdate);
            if  compare_matlab_versions(matlab_version, manifest_mat_version) == -1
                suffixes = [suffixes ['_compat.' manifest_mat_rdate]];
            end
        catch ME
            if ~strcmp(ME.identifier, 'Nirstorm:BadMatlabVersionString')
                rethrow(ME);
            end
        end
    end
end
if ~isempty(suffixes)
    s_extras = suffixes{1};
    for iextra=2:length(suffixes)
        s_extras = ['MANIFEST' s_extras ', ' suffixes{iextra}];
    end
    fprintf('Adding functions from %s to support matlab version %s.\n', ...
            s_extras, matlab_version);
end

end