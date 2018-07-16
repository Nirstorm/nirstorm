classdef WorkShopPerform2018Test < matlab.unittest.TestCase
    
    properties
        tmp_dir
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            tmpd = tempname;
            mkdir(tmpd);
            testCase.tmp_dir = tmpd;
            utest_bst_setup();
        end
    end
    
    methods(TestMethodTeardown)
        function tear_down(testCase)
            rmdir(testCase.tmp_dir, 's');
            utest_clean_bst();
        end
    end
    
    methods(Test)
        
        function test_downloader(testCase)
            global GlobalData;
            tmp_dir = testCase.tmp_dir;
            workshop_data_path = fullfile(nst_get_local_user_dir(), 'tutorials', 'Nirstorm_workshop_PERFORM_2018');
            if ~exist(workshop_data_path, 'dir')
                warning('Skipping test setup of workshop PERFORM 2018. Sample data directory not found in %s', ...
                        workshop_data_path);
                return; 
            end
            fluence_uncompress_dir = fullfile(workshop_data_path, 'fluences_for_OM');
            if exist(fluence_uncompress_dir , 'dir')
                warning('Uncompressed fluence folder exists, deleting it...');
                rmdir(fluence_uncompress_dir, 's');
            end
            
            if ~exist(workshop_data_path, 'dir')
                % TODO: check that all data files are there
                % data_fns =  process_nst_get_data_perform_2018('get_data_file_names');
                warning('Skip WorkshopPerform2018.test_downloader because data not found in %s', workshop_data_path); 
            end
            
            % Test behaviour when all is fine
            bst_dir = fullfile(tmpdir, 'bst_home');
            bst_process('CallProcess', ...
                        'process_nst_get_data_perform_2018', [], [], ...
                        'inputdir', {workshop_data_path, 'DIR'}, ...
                        'bst_dir', bst_dir);
                    
            installed_fns = [process_nst_get_data_perform_2018('get_colin27_installed_fn', bst_dir) ...
                             process_nst_get_data_perform_2018('get_OM_fluence_fns', bst_dir)];
                         
            fluence_uncompress_dir = fullfile(workshop_data_path, 'fluences_for_OM');
            testCase.assertTrue(exist(fluence_uncompress_dir, 'dir')==0);
            
            testCase.assertEmpty(GlobalData.lastestFullErrMsg);
            testCase.assertTrue(files_exist(installed_fns));
        end
        
    end
end

function flags = files_exist(fns)
flags = all(cellfun(@(fn) exist(fn, 'file')>0, fns));
end