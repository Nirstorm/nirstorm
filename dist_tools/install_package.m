function install_package(package_name, source_dir, target_dir, mode, extras, dry, target_matlab_ver)
% Install a package and manage version. From given package_name and a version 
% tag that must be defined in source_dir/VERSION, the function ensures 
% that only one package version is installed in the target directory.
% The set of files/folders to install are declared in the MANIFEST file that must
% be available in the current directory (see below "MANIFEST file format").
% Installation can copy files (windows/linux) or make symbolic links (linux only) 
% depending on the value of mode: 'copy' or 'link'. respectively.
% A version tag must be specified in source_dir/VERSION (see below 
% "VERSION file format")
% A uninstall script named uninstall_<package_name>.m is created. 
% It is used to clean a previous installation.
% If items to be installed already exist in target_dir disregarding
% uninstall script, the user is warned and they are backuped by prefixing
% "_backuped_by_<package_name>_<version_tag>_". The uninstallation script 
% will restore them.
%
% ## MANIFEST file format ##
% one path per line. Each path is relative to the folder of the MANIFEST file .
% Path can be a filename or a folder. If it's a folder then the whole 
% directory is copied or linked.
% Example:
%   my_func.m
%   script/my_script.m
%   data
%   
% ## VERSION file format ##
% contains only one string on the first line. Only alphanumerical
% characters and "." are allowed.
% IMPORTANT: version tag is case insentive, eg the version tage "stable_v1.2" 
% is considered the same as "Stable_V1.2".

% Args:
%     - package_name (str):
%         package name. Must be a valid matlab identifier. See function isvarname
%     - source_dir (str):
%         source directory which contains files to install.
%     - target_dir (str):
%         target installation folder.
%    [- mode (str):]
%         if 'copy': all files/folder are copied to target_dir
%         if 'link': symbolic links pointing to items in source_dir are
%                    created in target_dir.
%         Default is 'copy'.
%    [- extras (cell array of str):]
%         Tags to include extra the content of other MANIFEST files.
%         For each tag, a MANIFEST file named "MANIFEST.<tag>" must be
%         available in the source directory.
%    [- dry (boolean):]
%         If 1 then no file operations are displayed, not executed.
%         If 0 (default) then file operations are performed.
%    [- target_matlab_ver (str):]
%         Specify the target matlab version to consider for compatibility.
%         This is mainly for tests since it is always better to use
%         the current matlab version

if nargin < 4
    mode = 'copy';
end
if nargin < 5
    extras = {};
end
if nargin < 6
    dry = 0;
end

if nargin < 7
    target_matlab_ver = version(); 
end

check_inputs(package_name, source_dir, target_dir, mode, extras, dry);

uninstall_package(package_name, target_dir);
version_tag = get_version_tag(source_dir);
[install_operations, uninstall_operations] = resolve_file_operations(package_name, source_dir, target_dir, mode, extras, target_matlab_ver);
disp(['Installing ' package_name '--' version_tag ' to ' target_dir '...']);
execute_file_operations(install_operations, dry);
uninstall_script = fullfile(target_dir, ['uninstall_' package_name '.m']);
uninstall_header = sprintf('disp(''Uninstalling %s--%s from %s...'');', ...
                           package_name, version_tag, protect_path(target_dir));
make_file_operations_script(uninstall_operations, uninstall_script, uninstall_header, dry);
end

function make_file_operations_script(operations, script_fn, header, dry)

content = {header};
for iop=1:length(operations)
    operation = operations(iop);
    if ~isempty(operation.file1)
        if ~isfield(operation, 'dont_check_file1')
            content{end+1} = code_check_file_exists(operation.file1);
        else
            content{end+1} = sprintf('if exist(''%s'', ''file'')', operation.file1);
        end
        if ~isfield(operation, 'dont_check_file2') && ~isempty(operation.file2)
            content{end+1} = code_check_file_doesnt_exist(operation.file2);
        end
        switch operation.action
            case 'copy'
                content{end+1} = sprtinf('copyfile ''%s'' ''%s'';', operation.file1, operation.file2);
            case 'link'
                content{end+1} = code_check_not_windows();
                content{end+1} = sprintf('unix(''ln -s %s %s'')', operation.file1, operation.file2);
            case 'remove'              
                content{end+1} = strjoin_({sprintf('if ~isempty(strfind(computer, ''WIN'')) || unix(''test -L %s'')', operation.file1), ...
                                          sprintf('    if isdir(''%s'')', operation.file1), ...
                                          sprintf('        rmdir(''%s'', ''s'');', operation.file1), ...
                                                  '    else', ...
                                          sprintf('        delete(''%s'');', operation.file1), ...
                                                  '    end', ...
                                                  'else',...
                                          sprintf('    delete(''%s'');', operation.file1), ...       
                                                  'end',...
                                                  'try', ...
                                          sprintf('    rmdir(''%s'');', fileparts(operation.file1)),...
                                                  'catch', ...
                                                  'end'}, char(10));
            case 'move'
                content{end+1} = sprintf('movefile(''%s'', ''%s'');\n', operation.file1, operation.file2);
            otherwise
                throw(MException('DistPackage:BadOperation','Bad operation: %s', operation.action));
        end
        if isfield(operation, 'dont_check_file1')
            content{end+1} = 'end';
        end
    end
end
[rr, bfn, ext] = fileparts(script_fn);
func_header = sprintf('function %s()', bfn);
content = [func_header content 'end'];
content = sprintf('%s\n', strjoin_(content, char(10)));
if ~dry
    fout = fopen(script_fn, 'w');
    fprintf(fout, content);
    fclose(fout);
else
    fprintf(content);
end
end

function code = code_check_file_exists(fn)
code = {sprintf('if ~exist(fullfile(pwd, ''%s''), ''file'')', fn), ...
                '    dir_res = dir(fullfile(pwd));', ...
        sprintf('    if ~ismember(''%s'', {dir_res.name})', fn), ... 
        sprintf('        throw(MException(''DistPackage:FileNotFound'',''"%s" not found''));',fn),...
                '    else', ...
        sprintf('        warning(''DistPackage:BrokenLink'', [fullfile(pwd, ''%s'') ''seems to be a broken link'']);',fn),...        
                '    end',...
                'end'};
code = strjoin_(code, char(10));
end

function code = code_check_file_doesnt_exist(fn)
code = {sprintf('if exist(fullfile(pwd,''%s''), ''file'')', fn), ...
        sprintf('    throw(MException(''DistPackage:FileExists'',''File "%s" exists''));',fn),...
        'end'};
code = strjoin_(code, char(10));
end

function code = code_check_not_windows()
code = {'if ~isempty(strfind(computer, ''WIN''))', ...
        '    throw(MException(''DistPackage:BadOperation'', ''windows not supported''));', ...
        'end'};
code = strjoin_(code, char(10));
end

function [install_operations, uninstall_operations] = resolve_file_operations(package_name, source_dir, target_dir, mode, extras, target_matlab_ver)
compat_suffixes = [get_older_matlab_extras(source_dir, target_matlab_ver) ...
                   get_installable_extras(source_dir, target_matlab_ver)];
manifest_fns = [{fullfile(source_dir, 'MANIFEST')} ...
                cellfun(@(extra_tag) fullfile(source_dir, ['MANIFEST.' extra_tag]), ...
                        extras, 'UniformOutput', false) ...
                cellfun(@(suffix) fullfile(source_dir, ['MANIFEST' suffix]), ...
                        compat_suffixes, 'UniformOutput', false)];
backup_prefix = ['_backuped_by_' package_name '_' get_version_tag(source_dir) '_'];
iop = 1;
uop = 1;
for im=1:length(manifest_fns)
    manifest_fn = manifest_fns{im};
    if ~exist(manifest_fn, 'file')
       throw(MException('DistPackage:FileNotFound', [manifest_fn ' does not exist in root dir']));
    end
    source_rfns = read_manifest(manifest_fn);
    for ifn=1:length(source_rfns)
        source_fn = fullfile(source_dir, source_rfns{ifn});
        if ~exist(source_fn, 'file')
            throw(MException('DistPackage:FileNotFound', [protect_path(source_fn) ' does not exist in source dir']));
        end
        target_fn = fullfile(target_dir, source_rfns{ifn});
        if exist(target_fn, 'file')
            backup_rfn = add_fn_prefix(source_rfns{ifn}, backup_prefix);
            warning('DistPackage:ExistingTarget', ...
                    ['"' protect_path(source_rfns{ifn}) '" already exists in target directory. ' ...
                     'It will be backuped to "' backup_rfn '"']);
            install_operations(iop).file1 = target_fn;
            install_operations(iop).action = 'move';
            install_operations(iop).file2 = fullfile(target_dir, backup_rfn);
            iop = iop + 1;
        end
        % Safer to use absolute path when mode is link
        if java.io.File(source_fn).isAbsolute()
            abs_source_fn = source_fn;
        else
            abs_source_fn = fullfile(pwd, source_fn);
        end
        install_operations(iop).file1 = abs_source_fn;
        install_operations(iop).action = mode;
        install_operations(iop).file2 = target_fn;
        iop = iop + 1;

        uninstall_operations(uop).file1 = source_rfns{ifn};
        uninstall_operations(uop).action = 'remove';
        uninstall_operations(uop).file2 = '';
        uop = uop + 1;
        
        if exist(target_fn, 'file')
            uninstall_operations(uop).file1 = backup_rfn;
            uninstall_operations(uop).action = 'move';
            uninstall_operations(uop).file2 = source_rfns{ifn};
            uninstall_operations(uop).dont_check_file1 = 1;
            uop = uop + 1;
        end
    end
end
end

function fns = read_manifest(manifest_fn)
fns = cellfun(@(fn) strtrim(fn), strsplit_(fileread(manifest_fn), char(10)), 'UniformOutput', false);
fns = fns(~cellfun(@isempty, fns));
[root, ignore_bfn, ignore_ext] = fileparts(manifest_fn);
files_not_found = cellfun(@(fn) ~exist(fullfile(root, fn), 'file'), fns);
if any(files_not_found)
    fns_not_found = cellfun(@(fn) protect_path(fn), fns(files_not_found), 'UniformOutput', false);
    throw(MException('DistPackage:FileNotFound', ...
                     sprintf('Non-existing files from %s:\n%s', ...
                             protect_path(manifest_fn), strjoin_(fns_not_found, char(10)))));
end
end

function new_fn = add_fn_prefix(fn, prefix)
[rr, bfn, ext] = fileparts(fn);
new_fn = fullfile(rr, [prefix bfn  ext]);
end

function execute_file_operations(operations, dry)
for iop=1:length(operations)
    operation = operations(iop);
    if isdir(operation.file2) && exist(operation.file2, 'dir') || exist(operation.file2, 'file')
        throw(MException('DistPackage:TargetExists', ... 
                         ['Target ' protect_path(operation.file2) ' already exists. ' ...
                          'Installation aborted, consider manually cleaning target directory.']));
    end
    if ~isfield(operation, 'dont_check_file1') && (isdir(operation.file1) && ~exist(operation.file1, 'dir') || ~exist(operation.file1, 'file'))
        throw(MException('DistPackage:FileNotFound', ...
                         [protect_path(operation.file1) ' does not exist. ' ...
                          'Installation aborted, consider manually cleaning target directory.']));
    end
    if ~dry
        switch operation.action
            case 'copy'
                dest_folder = fileparts(operation.file2);
                if ~exist(dest_folder, 'dir')
                    mkdir(dest_folder);
                end
                copyfile(operation.file1, operation.file2);
            case 'link'
                dest_folder = fileparts(operation.file2);
                if ~exist(dest_folder, 'dir')
                    mkdir(dest_folder);
                end
                unix(['ln -s ' operation.file1 ' ' operation.file2]);
            case 'remove'
                if ~isempty(strfind(computer, 'WIN')) || unix(['test -L ' operation.file1])
                    if isdir(operation.file1)
                        rmdir(operation.file1, 's');
                    else
                        delete(operation.file1);
                    end
                else %symlink
                    delete(operation.file1);
                end
            case 'move'
                dest_folder = fileparts(operation.file2);
                if ~exist(dest_folder, 'dir')
                    mkdir(dest_folder);
                end
                movefile(operation.file1, operation.file2);
            otherwise
                throw(MException('DistPackage:BadOperation', ...
                                 ['Bad operation: ' operation.action]))
        end
    else
        switch operation.action
            case 'copy'
                disp(['copy ' protect_path(operation.file1) ' to ' protect_path(operation.file2)]);
            case 'link'
                disp(['link ' protect_path(operation.file1) ' to ' protect_path(operation.file2)]);
            case 'remove'
                disp(['remove ' protect_path(operation.file1)]);
            case 'move'
                disp(['move ' protect_path(operation.file1) ' to ' protect_path(operation.file2)]);
        end
    end
end

end

function version_tag = get_version_tag(root_dir)

version_fn = fullfile(root_dir, 'VERSION');
if ~exist(version_fn, 'file')
   throw(MException('DistPackage:FileNotFound', 'VERSION file does not exist in root dir'));
end

content = fileread(version_fn);
if ~isempty(content) && strcmp(sprintf(char(10)), content(end))
    content = content(1:(end-1));
end
version_tag = strtrim(content);
if isempty(regexp(version_tag, '^[a-z0-9_.]+$', 'once'))
    throw(MException('DistPackage:BadVersionTag', ...
                    ['Bad version tag in VERSION "' version_tag '". Must not be empty and only contain '...
                     'alphanumerical characters, dot or underscore']));
end
version_tag = lower(version_tag);
end

function check_inputs(package_name, source_dir, target_dir, mode, extras, dry)
if~isvarname(package_name)
    throw(MException('DistPackage:BadPackageName', ...
                     'Package name must be a valid matlab identifier'));
end

if ~exist(source_dir, 'dir')
   throw(MException('DistPackage:DirNotFound', ...
                    sprintf('source_dir "%s "does not exist', ...
                            protect_path(source_dir))));
end

if ~exist(target_dir, 'dir')
   mkdir(target_dir);
end 

if ~(strcmp(mode, 'copy') || strcmp(mode, 'link'))
    throw(MException('DistPackage:BadOption', 'mode can either be "copy" or "link"'));
end

if ~isempty(strfind(computer, 'WIN')) && strcmp(mode, 'link')
    throw(MException('DistPackage:BadOption', 'link mode only available for linux'));
end

if ~iscell(extras) || any(cellfun(@(e) ~ischar(e), extras))
    throw(MException('DistPackage:BadOption', 'extra must be a cell array of str'));
end

if ~(dry==0 || dry==1)
   throw(MException('DistPackage:BadOption', 'dry must be either 1 or 0'));
end

end

function toks = strsplit_(s, delimiter)
% For compatibility
try
    toks = strsplit(s, delimiter);
catch
    d = strtrim(delimiter);
    p = strfind(s, d);
    if ~isempty(p)                
        nt = numel(p) + 1;
        toks = cell(1, nt);
        sp = 1;
        dl = length(delimiter);
        for i=1:(nt-1)
            toks{i} = strtrim(s(sp:p(i)-1));
            sp = p(i) + dl;
        end         
        toks{nt} = strtrim(s(sp:end));
    else
        toks = {s};
    end        
end
end


function ppath = protect_path(path)
ppath = strrep(path, '\', '\\');
end
