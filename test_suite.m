function test_suite()
tic
test_install();
disp(['Successfully ran unit tests in ' num2str(toc) ' seconds.']);
end

function test_install()
test_package_install();
test_package_install_input_errors();
end

function test_package_install()
package_rel_files = {fullfile('script', 'script1.m'), ...
                     'func1.m', 'func2.m', fullfile('extra', 'func_extra.m'), ...
                     'dont_touch_me.sh', fullfile('please','dont_touch_me.m'),...
                     fullfile('dep', 'dep_func1.m'), ...
                     fullfile('dep', 'dep_func2.m')};
files_not_to_install = {'dont_touch_me.sh', fullfile('please','dont_touch_me.m')};
[package_files, root_src_dir] = create_tmp_files(package_rel_files);

% Create MANIFEST files
manifest_base_fn = fullfile(root_src_dir, 'MANIFEST');
base_rel_fns = {'script/script1.m', 'func1.m', 'func2.m', 'dep'};
write_file(manifest_base_fn, [base_rel_fns(1:2) '' '' base_rel_fns(3:end)]);

manifest_extra_fn = fullfile(root_src_dir, 'MANIFEST.extra');
extra_rel_fns = {'extra/func_extra.m'};
write_file(manifest_extra_fn, extra_rel_fns);

% Create VERSION file
version_fn = fullfile(root_src_dir, 'VERSION');
version_tag = '0.2.5';
write_file(version_fn, {version_tag});

% Create target directory
target_dir = tempname;
mkdir(target_dir);

% Test installation and uninstallation, copy mode
package_name = 'my_package';
uninstall_script_fn = fullfile(target_dir, ['uninstall_' package_name '.m']);
install_package(package_name, root_src_dir, target_dir, 'copy');
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(files_exist(target_dir, base_rel_fns)));
assert(all(~files_are_symlinks(target_dir, base_rel_fns)));
assert(all(~files_exist(target_dir, extra_rel_fns)));
assert(all(~files_exist(target_dir,files_not_to_install)));
assert(exist(uninstall_script_fn, 'file')>0);

uninstall_package(package_name, target_dir);
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(~files_exist(target_dir, base_rel_fns)));
assert(~exist(uninstall_script_fn, 'file')>0);

% Test installation and uninstallation, link mode
install_package(package_name, root_src_dir, target_dir, 'link');
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(files_exist(target_dir, base_rel_fns)));
assert(all(files_are_symlinks(target_dir, base_rel_fns)));
assert(all(~files_exist(target_dir, extra_rel_fns)));
assert(all(~files_exist(target_dir,files_not_to_install)));
assert(exist(uninstall_script_fn, 'file')>0);

uninstall_package(package_name, target_dir);
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(~files_exist(target_dir, base_rel_fns)));
assert(~exist(uninstall_script_fn, 'file')>0);

% Test installation with already existing items in target_dir -> should be
% backuped
existing_target_fn = fullfile(target_dir, 'func1.m');
existing_target_backup_fn = add_fn_prefix(existing_target_fn, ...
                                          ['_backuped_by_' package_name '_' version_tag '_']);
fout = fopen(existing_target_fn, 'w');
fclose(fout);
install_package(package_name, root_src_dir, target_dir, 'link');
assert(exist(existing_target_backup_fn, 'file')>0);
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(files_exist(target_dir, base_rel_fns)));
assert(exist(uninstall_script_fn, 'file')>0);

uninstall_package(package_name, target_dir);
assert(exist(existing_target_fn, 'file')>0);
assert(exist(existing_target_backup_fn, 'file')==0);
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(~files_exist(target_dir, [base_rel_fns(1) base_rel_fns(3:end)])));
assert(~exist(uninstall_script_fn, 'file')>0);

rmdir(target_dir, 's');
rmdir(root_src_dir, 's');
end

function test_package_install_input_errors()
package_rel_files = {fullfile('script', 'script1.m'), ...
                     'func1.m', 'func2.m', fullfile('extra', 'func_extra.m'), ...
                     'dont_touch_me.sh', fullfile('please','dont_touch_me.m'),...
                     fullfile('dep', 'dep_func1.m'), ...
                     fullfile('dep', 'dep_func2.m')};
[package_files, root_src_dir] = create_tmp_files(package_rel_files);

% Create MANIFEST files
manifest_base_fn = fullfile(root_src_dir, 'MANIFEST');
base_rel_fns = {'script/script1.m', 'func1.m', 'func2.m', 'dep'};
write_file(manifest_base_fn, base_rel_fns);

manifest_bad_fn = fullfile(root_src_dir, 'MANIFEST.bad');
write_file(manifest_bad_fn, {'func_crowe.m', 'dir_crowe'});

% Create VERSION file
version_fn = fullfile(root_src_dir, 'VERSION');
version_tag = '0.2.5';
write_file(version_fn, {version_tag});

% Create target directory
target_dir = tempname;
mkdir(target_dir);

% Test installation with bad MANIFEST
package_name = 'my_package';
uninstall_script_fn = fullfile(target_dir, ['uninstall_' package_name '.m']);

exception_caught = 0;
try
    install_package(package_name, root_src_dir, target_dir, 'copy', {'bad'});
catch ME
    assert(strcmp(ME.identifier, 'DistPackage:FileNotFound'))
    assert(~isempty(strfind(ME.message, 'func_crowe.m')));
    assert(~isempty(strfind(ME.message, 'dir_crowe')));
    exception_caught = 1;
end
assert(all(files_exist(root_src_dir, package_rel_files)));
assert(all(~files_exist(target_dir, base_rel_fns)));
assert(~exist(uninstall_script_fn, 'file')>0);
assert(exception_caught==1);
end

function flags = files_exist(root_dir, rel_fns)
flags = cellfun(@(rfn) exist(fullfile(root_dir, rfn), 'file')>0, rel_fns);
end

function flags = files_are_symlinks(root_dir, rel_fns)
if ~isempty(strfind(computer, 'WIN'))
    flags = zeros(1, length(rel_fns))==0;
else
    flags = cellfun(@(rfn) unix(['test -L ' fullfile(root_dir, rfn)])==0, rel_fns); 
end
end

function write_file(fn, lines)
folder = fileparts(fn);
if ~exist(folder, 'dir')
    mkdir(folder);
end
fout = fopen(fn, 'w');
fprintf(fout, strjoin(lines, '\n'));
fclose(fout);
end

function [tmp_fns, tmp_dir] = create_tmp_files(rel_fns)
% Create empty temporary files from the given list of relative pathes

tmp_dir = tempname;
mkdir(tmp_dir);
tmp_fns = cell(1, length(rel_fns));
for ifn=1:length(rel_fns)
    tmp_fn = fullfile(tmp_dir, rel_fns{ifn});
    dest_folder = fileparts(tmp_fn);
    if ~exist(dest_folder, 'dir')
        mkdir(dest_folder);
    end
    fout = fopen(tmp_fn, 'w');
    fclose(fout);
    tmp_fns{ifn} = tmp_fn;
end

end

function new_fn = add_fn_prefix(fn, prefix)
[rr, bfn, ext] = fileparts(fn);
new_fn = fullfile(rr, [prefix bfn ext]);
end
