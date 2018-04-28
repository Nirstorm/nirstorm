classdef InstallSourceTest < matlab.unittest.TestCase
    methods(Test)
        function test_matlab_version_comparison(testCase)
            exception_caught = 0;
            try
                compare_matlab_versions('8.1', 9.1);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, ...
                                    'Nirstorm:compare_matlab_versions:invalidInput'));
                 exception_caught = 1;
            end
            testCase.assertTrue(exception_caught==1);
            
            testCase.assertEqual(compare_matlab_versions('8.1', '9.1.1'), -1);
            testCase.assertEqual(compare_matlab_versions('8.1', '9.1'), -1);
            testCase.assertEqual(compare_matlab_versions('8.1', '9'), -1);
            testCase.assertEqual(compare_matlab_versions('7.11.1', '7.2'), 1);
            testCase.assertEqual(compare_matlab_versions('7.11.1', '7.11.1'), 0);
        end
        function test_package_install(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture

            addpath('../dist_tools');
            
            testCase.applyFixture(SuppressedWarningsFixture('DistPackage:FileNotFound'));
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
            install_package(package_name, root_src_dir, target_dir, 'copy', {}, 0);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(all(~files_are_symlinks(target_dir, base_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir, extra_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir,files_not_to_install)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            uninstall_package(package_name, target_dir);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);

            % Test installation and uninstallation, link mode
            install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(all(files_are_symlinks(target_dir, base_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir, extra_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir,files_not_to_install)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            uninstall_package(package_name, target_dir);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);

            % Test installation with already existing items in target_dir 
            % -> should be backuped
            existing_target_fn = fullfile(target_dir, 'func1.m');
            existing_target_backup_fn = add_fn_prefix(existing_target_fn, ...
                                                      ['_backuped_by_' package_name '_' version_tag '_']);
            fout = fopen(existing_target_fn, 'w');
            fprintf(fout, 'blah');
            fclose(fout);
            install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            uninstall_package(package_name, target_dir);
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')==0);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, [base_rel_fns(1) base_rel_fns(3:end)])));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);

            % Test installation with backups items in target_dir 
            % that got deleted. Uninstall should not complain 
            existing_target_fn = fullfile(target_dir, 'func1.m');
            existing_target_backup_fn = add_fn_prefix(existing_target_fn, ...
                                                      ['_backuped_by_' package_name '_' version_tag '_']);
            fout = fopen(existing_target_fn, 'w');
            fprintf(fout, 'blah');
            fclose(fout);
            install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            delete(existing_target_backup_fn);
            uninstall_package(package_name, target_dir);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            
            rmdir(target_dir, 's');
            rmdir(root_src_dir, 's');
        end
        
        function test_package_install_errors(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture

            addpath('../dist_tools');
            
            testCase.applyFixture(SuppressedWarningsFixture('DistPackage:FileNotFound'));
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
            version_tag = sprintf('0.2.5');
            write_file(version_fn, {version_tag});

            % Create target directory
            target_dir = tempname;
            mkdir(target_dir);

            % Test installation with bad MANIFEST
            package_name = 'my_package';
            uninstall_script_fn = fullfile(target_dir, ['uninstall_' package_name '.m']);

            exception_caught = 0;
            try
                install_package(package_name, root_src_dir, target_dir, 'copy', {'bad'}, 0);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'DistPackage:FileNotFound'))
                testCase.assertTrue(~isempty(strfind(ME.message, 'func_crowe.m')));
                testCase.assertTrue(~isempty(strfind(ME.message, 'dir_crowe')));
                exception_caught = 1;
            end
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            testCase.assertTrue(exception_caught==1);

            % Test with already existing backuped file
            existing_target_fn = fullfile(target_dir, 'func1.m');
            existing_target_backup_fn = add_fn_prefix(existing_target_fn, ...
                                                      ['_backuped_by_' package_name '_' version_tag '_']);
            fout = fopen(existing_target_backup_fn, 'w');
            fprintf(fout, 'blah');
            fclose(fout);
            fout = fopen(existing_target_fn, 'w');
            fprintf(fout, 'meh');
            fclose(fout);
            exception_caught = 0;
            try
                install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'DistPackage:TargetExists'))
                testCase.assertTrue(~isempty(strfind(ME.message, existing_target_backup_fn)));
                exception_caught = 1;
            end
            testCase.assertTrue(exception_caught==1);
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            
            % No version file
            delete(version_fn)
            exception_caught = 0;
            try
                install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'DistPackage:FileNotFound'))
                testCase.assertTrue(~isempty(strfind(ME.message, 'VERSION file')));
                exception_caught = 1;
            end
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            testCase.assertTrue(exception_caught==1);
            
            % Bad version tags
            version_tag = '|@m reA11y b@@@@@d';
            write_file(version_fn, {version_tag});
            exception_caught = 0;
            try
                install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'DistPackage:BadVersionTag'))
                testCase.assertTrue(~isempty(strfind(ME.message, version_tag)));
                exception_caught = 1;
            end
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            testCase.assertTrue(exception_caught==1);
            
            version_tag = '';
            write_file(version_fn, {version_tag});
            exception_caught = 0;
            try
                install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'DistPackage:BadVersionTag'))
                exception_caught = 1;
            end
            testCase.assertTrue(all(files_exist(root_src_dir, package_rel_files)));
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            testCase.assertTrue(exception_caught==1);
            
            rmdir(target_dir, 's');
            rmdir(root_src_dir, 's');
        end
        
        function test_uninstall_error(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture

            addpath('../dist_tools');
            
            testCase.applyFixture(SuppressedWarningsFixture('DistPackage:FileNotFound'));
            package_rel_files = {fullfile('script', 'script1.m'), ...
                     'func1.m', 'func2.m', fullfile('extra', 'func_extra.m'), ...
                     'dont_touch_me.sh', fullfile('please','dont_touch_me.m'),...
                     fullfile('dep', 'dep_func1.m'), ...
                     fullfile('dep', 'dep_func2.m')};
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

            % Test uninstallation where installed file has been deleted
            package_name = 'my_package';
            existing_target_fn = fullfile(target_dir, 'func1.m');
            existing_target_backup_fn = add_fn_prefix(existing_target_fn, ...
                                                      ['_backuped_by_' package_name '_' version_tag '_']);
            fout = fopen(existing_target_backup_fn, 'w');
            fprintf(fout, 'blah');
            fclose(fout);
            
            install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            delete(existing_target_backup_fn);
            delete(existing_target_fn);
            exception_caught = 0;
            try
                uninstall_package(package_name, target_dir);
            catch ME
                testCase.assertTrue(strcmp(ME.identifier, 'DistPackage:FileNotFound'))
                testCase.assertTrue(~isempty(strfind(ME.message, 'func1.m')));
                exception_caught = 1;
            end
            testCase.assertTrue(exception_caught==1);
            
            rmdir(target_dir, 's');
            rmdir(root_src_dir, 's');
        end
    end
end

function new_fn = add_fn_prefix(fn, prefix)
[rr, bfn, ext] = fileparts(fn);
new_fn = fullfile(rr, [prefix bfn ext]);
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

function write_file(fn, lines)
folder = fileparts(fn);
if ~exist(folder, 'dir')
    mkdir(folder);
end
fout = fopen(fn, 'w');
fprintf(fout, strjoin(lines, '\n'));
fclose(fout);
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
