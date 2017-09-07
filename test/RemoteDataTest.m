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
        
    end
end

function flag = nst_url_equal(u1, u2)
flag = strcmp(normalize_url(u1), normalize_url(u2));
end

function u = normalize_url(u)
u = regexprep(u,'//+', '/');
end
