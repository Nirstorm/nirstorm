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
        
        function test_run_proc_single_output(testCase)
            global GlobalData
            
            % Prepare data
            y = [1:10 ; 101:110];
            time = (0:(size(y,2)-1))*0.1;
            
            bst_create_test_subject();
            
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            sFilesOut = nst_run_bst_proc('dummy_resampled', 0, 'process_resample', sRaw, [], 'freq', 5);
            
            % TODO: test that output has correct Comment
            % TODO: test that resampling was applied
            
            %% Call proc again
            sFilesOut = nst_run_bst_proc('dummy_resampled', 0, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertEmpty(sFilesOut);
            % TODO: test that output 'dummy_resampled' still exists
            
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, ...
                                   'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,4}, ...
                                   'nst_run_bst_proc: skipped (output(s) found)');                    
        end
        
        function test_run_proc_multiple_outputs(testCase)
            
        end

        function test_run_proc_dont_redo(testCase)
            
        end

        function test_run_proc_force_redo(testCase)
            
        end
        
    end
end