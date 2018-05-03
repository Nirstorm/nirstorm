classdef InstallSourceTest < matlab.unittest.TestCase
    properties
        tmp_dir
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
        end
    end
    
    methods(TestMethodTeardown)
        function tear_down(testCase)
            rmdir(testCase.tmp_dir, 's');
        end
    end
    
    methods(Test)
                
        function test_older_matlab_install(testCase)
            create_tmp_files(testCase, {'MANIFEST_compat.R2008a', 'MANIFEST_compat.R2009b', 'MANIFEST_compat.R2014b'});
            
            suffixes = get_older_matlab_extras(testCase.tmp_dir, rdate_to_version('R2008a'));
            testCase.assertTrue(all(ismember(suffixes, {'_compat.R2009b', '_compat.R2014b'})) && length(suffixes)==2);
            
            suffixes = get_older_matlab_extras(testCase.tmp_dir, rdate_to_version('R2007a'));
            testCase.assertTrue(all(ismember(suffixes, {'_compat.R2008a', '_compat.R2009b', '_compat.R2014b'})) && length(suffixes)==3);
            
            suffixes = get_older_matlab_extras(testCase.tmp_dir, rdate_to_version('R2015a'));
            testCase.assertTrue(isempty(suffixes));
            
%             To test extras on the current matlab version (number of extras may vary):         
%             extras = get_older_matlab_extras(testCase.tmp_dir);
%             testCase.assertTrue(isempty(extras));

        end
        
        function test_specific_matlab_install(testCase)
            manifest_files = create_tmp_files(testCase, {'MANIFEST_only_from.R2008a', ...
                                                         'MANIFEST_only_from.R2009b', ...
                                                         'MANIFEST_only_from.R2014b'});
            write_file(manifest_files{1}, {'func1_only_from_R2008a.m', 'func2_only_from_R2008a.m'});
            write_file(manifest_files{2}, {'func1_only_from_R2009b.m', 'func2_only_from_R2009b.m'});
            write_file(manifest_files{3}, {'func1_only_from_R2014b.m', 'func2_only_from_R2014b.m'});
            
            suffixes = get_installable_extras(testCase.tmp_dir, rdate_to_version('R2008a'));
            testCase.assertTrue(all(ismember(suffixes, {'_only_from.R2008a'})) && length(suffixes)==1);
            
            suffixes = get_installable_extras(testCase.tmp_dir, rdate_to_version('R2007a'));
            testCase.assertTrue(isempty(suffixes));
            
            suffixes = get_installable_extras(testCase.tmp_dir, rdate_to_version('R2015a'));
            testCase.assertTrue(all(ismember(suffixes, {'_only_from.R2008a', '_only_from.R2009b', '_only_from.R2014b'})) && length(suffixes)==3);            
        end
        
        
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
        
        
        function test_matlab_version_specifics(testCase)
            base_rel_fns = { 'func1.m', 'func_R2015b_alternate.m', ...
                             'func_only_for_R2016a.m'};
            package_rel_files = cellfun(@(rfn) fullfile('src_root', rfn), base_rel_fns, 'UniformOutput', false);
            create_tmp_files(testCase, package_rel_files);
            root_src_dir = fullfile(testCase.tmp_dir, 'src_root');

            % Create MANIFEST files
            manifest_base_fn = fullfile(root_src_dir, 'MANIFEST');
            write_file(manifest_base_fn, {'func1.m'});

            manifest_compat_fn = fullfile(root_src_dir, 'MANIFEST_compat.R2015b');
            write_file(manifest_compat_fn, {'func_R2015b_alternate.m'});

            manifest_vspec_fn = fullfile(root_src_dir, 'MANIFEST_only_from.R2016a');
            write_file(manifest_vspec_fn, {'func_only_for_R2016a.m'});
            
            % Create VERSION file
            version_fn = fullfile(root_src_dir, 'VERSION');
            version_tag = '0.2.5';
            write_file(version_fn, {version_tag});

            % Create target directory
            target_dir = fullfile(testCase.tmp_dir, 'target');
            mkdir(target_dir);

            % Test installation and uninstallation, copy mode
            package_name = 'my_package';
            matlab_version = '7.10'; % R0210a
            install_package(package_name, root_src_dir, target_dir, 'copy', {}, 0, matlab_version);
            testCase.assertTrue(all(files_exist(target_dir, {base_rel_fns{1:2}})));
            testCase.assertTrue(all(~files_exist(target_dir, {base_rel_fns{3}})));
            uninstall_package(package_name, target_dir);
            
            matlab_version = '9.4'; % R0218b
            install_package(package_name, root_src_dir, target_dir, 'copy', {}, 0, matlab_version);
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns([1 3]))));
            testCase.assertTrue(all(~files_exist(target_dir, {base_rel_fns{2}})));
            uninstall_package(package_name, target_dir);
        end
        
        function test_package_install(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture

            addpath('../dist_tools');
            
            testCase.applyFixture(SuppressedWarningsFixture('DistPackage:FileNotFound'));
                        
            package_rel_files = {fullfile('src_root', 'script', 'script1.m'), ...
                                 fullfile('src_root', 'func1.m'), ...
                                 fullfile('src_root', 'func2.m'), ...
                                 fullfile('src_root', 'extra', 'func_extra.m'), ...
                                 fullfile('src_root', 'dont_touch_me.sh'),...
                                 fullfile('src_root', 'please','dont_touch_me.m'),...
                                 fullfile('src_root', 'dep', 'dep_func1.m'), ...
                                 fullfile('src_root', 'dep', 'dep_func2.m')};
            files_not_to_install = {fullfile('src_root', 'dont_touch_me.sh'), ...
                                    fullfile('src_root', 'please','dont_touch_me.m')};
            create_tmp_files(testCase, package_rel_files);
            root_src_dir = fullfile(testCase.tmp_dir, 'src_root');

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
            target_dir = fullfile(testCase.tmp_dir, 'target');
            mkdir(target_dir);

            % Test installation and uninstallation, copy mode
            package_name = 'my_package';
            uninstall_script_fn = fullfile(target_dir, ['uninstall_' package_name '.m']);
            install_package(package_name, root_src_dir, target_dir, 'copy', {}, 0);
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(all(~files_are_symlinks(target_dir, base_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir, extra_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir,files_not_to_install)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            uninstall_package(package_name, target_dir);
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);

            % Test installation and uninstallation, link mode
            install_package(package_name, root_src_dir, target_dir, 'link', {}, 0);
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(all(files_are_symlinks(target_dir, base_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir, extra_rel_fns)));
            testCase.assertTrue(all(~files_exist(target_dir,files_not_to_install)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            uninstall_package(package_name, target_dir);
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            uninstall_package(package_name, target_dir);
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')==0);
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(all(files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(exist(uninstall_script_fn, 'file')>0);

            delete(existing_target_backup_fn);
            uninstall_package(package_name, target_dir);
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(all(~files_exist(target_dir, base_rel_fns)));
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            
        end
        
        function test_package_install_errors(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture

            addpath('../dist_tools');
            
            testCase.applyFixture(SuppressedWarningsFixture('DistPackage:FileNotFound'));
            
            package_rel_files = {fullfile('src_root', 'script', 'script1.m'), ...
                                 fullfile('src_root', 'func1.m'), ...
                                 fullfile('src_root', 'func2.m'), ...
                                 fullfile('src_root', 'extra', 'func_extra.m'), ...
                                 fullfile('src_root', 'dont_touch_me.sh'),...
                                 fullfile('src_root', 'please','dont_touch_me.m'),...
                                 fullfile('src_root', 'dep', 'dep_func1.m'), ...
                                 fullfile('src_root', 'dep', 'dep_func2.m')};            
            create_tmp_files(testCase, package_rel_files);
            root_src_dir = fullfile(testCase.tmp_dir, 'src_root');

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
            target_dir = fullfile(testCase.tmp_dir, 'target');
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
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
            testCase.assertTrue(all(files_exist(testCase.tmp_dir, package_rel_files)));
            testCase.assertTrue(exist(existing_target_fn, 'file')>0);
            testCase.assertTrue(exist(existing_target_backup_fn, 'file')>0);
            testCase.assertTrue(~exist(uninstall_script_fn, 'file')>0);
            testCase.assertTrue(exception_caught==1);
            
        end
        
        function test_uninstall_error(testCase)
            import matlab.unittest.fixtures.SuppressedWarningsFixture

            addpath('../dist_tools');
            
            testCase.applyFixture(SuppressedWarningsFixture('DistPackage:FileNotFound'));
            package_rel_files = {fullfile('src_root', 'script', 'script1.m'), ...
                                 fullfile('src_root', 'func1.m'), ...
                                 fullfile('src_root', 'func2.m'), ...
                                 fullfile('src_root', 'extra', 'func_extra.m'), ...
                                 fullfile('src_root', 'dont_touch_me.sh'),...
                                 fullfile('src_root', 'please','dont_touch_me.m'),...
                                 fullfile('src_root', 'dep', 'dep_func1.m'), ...
                                 fullfile('src_root', 'dep', 'dep_func2.m')};            
            create_tmp_files(testCase, package_rel_files);
            root_src_dir = fullfile(testCase.tmp_dir, 'src_root');
            
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
            target_dir = fullfile(testCase.tmp_dir, 'target');
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
        end
    end
end

function new_fn = add_fn_prefix(fn, prefix)
[rr, bfn, ext] = fileparts(fn);
new_fn = fullfile(rr, [prefix bfn ext]);
end

function [tmp_fns] = create_tmp_files(testCase, rel_fns)
% Create empty temporary files from the given list of relative pathes

tmp_fns = cell(1, length(rel_fns));
for ifn=1:length(rel_fns)
    tmp_fn = fullfile(testCase.tmp_dir, rel_fns{ifn});
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
