function extras = get_older_matlab_extras(manifest_dir, matlab_version)
% Append extra installation scenarios for older versions of matlab.
% Add extras corresponding to MANIFEST.<mat_release> where <release_date>
% is more recent than given matlab.
%
% Args:
%    - manifest_dir (str):
%         Directory where to search for MANIFEST.*
%    [- matlab_version (str)]:
%         matlab version string that must be supported.
%         Default is current matlab version
%
% Outputs:
%     - extras (cell of str):
%         List of extra scenario labels corresponding to MANIFEST.extra
%         files that will install functions to support an older matlab
%

if nargin < 2
    matlab_version = version();
end

items = dir(manifest_dir);
extras = {};
for iitem=1:length(items)
    item = items(iitem);
    [root, basename, ext] = fileparts(item.name);
    try
        manifest_mat_rdate = ext(2:end); % discard dot
        manifest_mat_version = rdate_to_version(manifest_mat_rdate); 
        if ~item.isdir && strcmp(basename, 'MANIFEST') && ...
                compare_matlab_versions(matlab_version, manifest_mat_version) == -1
            extras = [extras manifest_mat_rdate];
        end
    catch ME
        if ~strcmp(ME.identifier, 'Nirstorm:BadMatlabVersionString')
            rethrow(ME);
        end
    end
end
if ~isempty(extras)
    s_extras = extras{1};
    for iextra=2:length(extras)
        s_extras = [s_extras ', ' extras{iextra}];
    end
    fprintf('Adding functions for matlab versions %s to support current version %s.\n', ...
            s_extras, matlab_version);
end

end