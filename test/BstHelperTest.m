classdef BstHelperTest < matlab.unittest.TestCase
    
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
        
        function test_run_proc_output_mismatch(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            sFilesOut = nst_run_bst_proc({'dummy_resampled', 'dummy_dummy'}, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertTrue(~isempty(sFilesOut));
            testCase.assertTrue(~isempty(strfind(GlobalData.lastestFullErrMsg, ...
                                'Expected 2 outputs but process produced 1')));                  
        end
        
        function test_run_proc_duplicate_outputs(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            output_name = 'dummy_resampled';
            sFilesOut = bst_process('CallProcess', 'process_resample', sRaw, [], 'freq', 5);
            bst_process('CallProcess', 'process_set_comment', sFilesOut, [], ...
                        'tag', output_name, 'isindex', 0);
            
            sFilesOut = bst_process('CallProcess', 'process_resample', sRaw, [], 'freq', 2);
            bst_process('CallProcess', 'process_set_comment', sFilesOut, [], ...
                        'tag', output_name, 'isindex', 0);
                    
                    
            sFilesOut = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);
            testCase.assertEmpty(sFilesOut);
            testCase.assertTrue(~isempty(strfind(GlobalData.lastestFullErrMsg, ...
                                'Cannot safely manage unique outputs. Found duplicate items: dummy_resampled')));

        end
        
        
        function test_run_proc_dont_redo(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            output_name = 'dummy_resampled';
            sFilesOut = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');

            testCase.assertMatches(sFilesOut.Comment, output_name);
            resampled_data = in_bst_data(sFilesOut.FileName);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            sFilesOut = nst_run_bst_proc('dummy_resampled', 0, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Set comment');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,4}, ...
                                   'Skipped execution of process_resample. Outputs found.');                    
        end

        function test_run_proc_force_redo(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            output_name = 'dummy_resampled';
            sFilesOut = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');

            testCase.assertMatches(sFilesOut.Comment, output_name);
            resampled_data = in_bst_data(sFilesOut.FileName);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            sFilesOut = nst_run_bst_proc('dummy_resampled', 1, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut.FileName), 'file')==2);
            
            testCase.assertTrue(~isempty(strfind(GlobalData.ProcessReports.Reports{end-3,4}, ...
                                                 'Force redo - removing previous result(s):')));
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-2,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-2,2}.Comment, 'Delete files');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');

        end
        
        
        
        function test_run_proc_multiple_outputs(testCase)
            
        end
        
    end
end