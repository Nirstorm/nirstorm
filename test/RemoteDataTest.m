classdef RemoteDataTest < matlab.unittest.TestCase
    
    properties
        tmp_dir
    end
    
    methods(TestMethodSetup)
        function create_tmp_dir(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
        end
    end
    
    methods(TestMethodTeardown)
        function rm_tmp_dir(testCase)
            rmdir(testCase.tmp_dir, 's');
        end
    end
    
    methods(Test)
        function test_file_request(testCase)
            %% Test download only when file not available locally
            % Prepare local files
            utest_dir = fullfile(nst_get_local_user_dir(), 'unittest');
            if ~exist(utest_dir, 'dir')
                mkdir(utest_dir);
            end
            fout = fopen(fullfile(nst_get_local_user_dir(), 'unittest', 'dummy.txt'), 'w');
            fclose(fout);
            data_fns = {{'unittest', 'dummy.txt'}, {'unittest', 'motor_data', 'fiducials.txt'}, ...
                        {'unittest', 'motor_data', 'optodes.txt'}};
            for ifns=2:length(data_fns) % make sure some files are not locally available
                fn = fullfile(nst_get_local_user_dir(), strjoin(data_fns{ifns}, filesep));
                if exist(fn, 'file')
                    delete(fn);
                end
            end
            % Request files
            [local_fns, downloaded_fns] = nst_request_files(data_fns, 0);
            
            % Check results
            assert(all(cellfun(@(fn) ~isempty(strfind(fn, nst_get_local_user_dir())), local_fns)));
            assert(all(cellfun(@(fn) exist(fn, 'file'), local_fns)));
            assert(length(downloaded_fns) == 2);
            assert(strcmp(downloaded_fns{1}, strjoin({nst_get_repository_url(), strjoin(data_fns{2}, '/')}, '/')));
            assert(strcmp(downloaded_fns{2}, strjoin({nst_get_repository_url(), strjoin(data_fns{3}, '/')}, '/')));
            
            %% Test error when some requested file not available locally nor remotely
            data_fns = {{'unittest', 'dummy.txt'}, {'unittest', 'malcolm_crowe'}} ;
            try
                nst_request_files(data_fns, 0);
            catch ME
                testCase.assertMatches(ME.identifier, 'NIRSTORM:RemoteFilesNotFound');
                files_not_found = strsplit(ME.message, '\n');
                files_not_found = files_not_found(2:end);
                assert(nst_url_equal(files_not_found{1}, ...
                    strjoin({nst_get_repository_url(), strjoin(data_fns{2}, '/')}, '/')));
            end
            
        end
        
        function test_repository_sanity(testCase)
            errors = check_remote_repo_sanity_rec(nst_get_repository_url(), '', {});
            disp('Checking consistency between local and remote files (may take a bit of time).');
            warnings = check_local_repo_sanity_rec(nst_get_local_user_dir(), nst_get_repository_url(), {}, {});
            disp('Done checking consistency between local and remote files.');
            
            if ~isempty(warnings)
                warning(strjoin(warnings, '\n'));
            end
            
            if ~isempty(errors)
                error(strjoin(errors, '\n'));
            end
        end
        
    end
end



function errors = check_local_repo_sanity_rec(local_root, url_root, subdirs, errors)

% disp('checking: ');
% disp(subdirs);
subdir = strjoin(subdirs, filesep);
local_items = dir(fullfile(local_root, subdir));
for ii=1:length(local_items)
    if ~local_items(ii).isdir
        disp(['checking ' local_items(ii).name]);
        local_file = fullfile(local_root, subdir, local_items(ii).name);
        remote_file = strjoin({url_root, strjoin(subdirs, '/'), local_items(ii).name}, '/');
        if ~nst_url_exists(remote_file)
            errors = [errors ['Local file "' local_file ...
                              '" not remotely available at "' ...
                              remote_file '".']];
        end
    elseif ~strcmp(local_items(ii).name, '.') && ~strcmp(local_items(ii).name, '..')
        errors = check_local_repo_sanity_rec(local_root, url_root, ...
                                             [subdirs local_items(ii).name], ...
                                             errors);
    end
end

end

function errors = check_remote_repo_sanity_rec(url_root, subdir, errors)
listing_tmp_fn = [tempname '.html'];
if strcmp(subdir, './') || strcmp(subdir, '../')
    return;
end
websave(listing_tmp_fn, [url_root subdir]);
listing_html = fileread(listing_tmp_fn);
delete(listing_tmp_fn);

files = regexp(listing_html, 'B  <a href="(.*?)">', 'tokens');
if ~isempty(files)
    listing_file_sizes_url = [url_root subdir 'file_sizes.csv'];
    if ~nst_url_exists(listing_file_sizes_url)
        errors = [errors ['Missing file size listing in ' [url_root subdir]]];
    else
        listing_file_sizes_fn = [tempname '.csv'];
        websave(listing_file_sizes_fn, listing_file_sizes_url);
        file_sizes = readtable(listing_file_sizes_fn, 'Delimiter', ',');
        for ifn=1:length(files)
            cur_file = files{ifn}{1};
            if ~strcmp(cur_file, 'file_sizes.csv') && ~ismember(cur_file, file_sizes.file_name)
                errors = [errors ['Missing file size of ' [subdir cur_file]]];
            end
        end
    end
end
subdirs = regexp(listing_html, '&lt;DIR&gt;  <a href="(.*?)">', 'tokens');
for isubdir=1:length(subdirs)
    errors = check_remote_repo_sanity_rec([url_root subdir], subdirs{isubdir}{1}, errors);
end

end

function flag = nst_url_exists(url)
jurl = java.net.URL(url);
conn = openConnection(jurl);
status = getResponseCode(conn);
flag = status ~= 404;
end


function flag = nst_url_equal(u1, u2)
flag = strcmp(normalize_url(u1), normalize_url(u2));
end

function u = normalize_url(u)
u = regexprep(u,'//+', '/');
end
