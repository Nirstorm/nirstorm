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
            
            resampled_data = in_bst_data(sFilesOut);
            testCase.assertMatches(resampled_data.Comment, output_name);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            sFilesOut = nst_run_bst_proc('dummy_resampled', 1, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-3,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-3,2}.Comment, 'Delete files');
            testCase.assertTrue(~isempty(strfind(GlobalData.ProcessReports.Reports{end-2,4}, ...
                                                 'Force redo - removed previous result(s):')));
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Set comment');

        end

        
        function test_run_proc_new_cond_dont_redo(testCase)
            global GlobalData
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);
                                    
            %% 1st proc call
            output_cond = 'new_cond';
            output_comment = 'dummy_resampled';
            output_name = [output_cond '/' output_comment];
            sFilesOut = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);

            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-2,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-2,2}.Comment, 'Resample');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Set comment');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Move files');
            prev_report_nb_items = length(GlobalData.ProcessReports.Reports);
            
            testCase.assertNotEmpty(nst_get_bst_func_files('test_subject', output_cond, output_comment));
            
            resampled_data = in_bst_data(sFilesOut);
            testCase.assertMatches(resampled_data.Comment, output_comment);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            sFilesOut = nst_run_bst_proc(output_name, 0, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            testCase.assertEqual(length(GlobalData.ProcessReports.Reports), prev_report_nb_items+1);
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Move files'); %from previous call
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,4}, ...
                                   'Skipped execution of process_resample. Outputs found.');                    
        end 
        
        
        function test_head_model_force_redo(testCase)
            %TODO
        end
        
        function test_head_model_redo(testCase)
            global GlobalData;
            
            %% Prepare data
            y = [1:10 ; 101:110];
            dt = 0.1;
            time = (0:(size(y,2)-1)) * dt;
            
            bst_create_test_subject();
            
            sRaw = bst_create_nirs_data('dummy_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]); 
                                    
            %% 1st proc call
            output_name = 'head_model';
            sFilesOut = nst_run_bst_proc(output_name, 0, 'process_nst_import_head_model', ...
                                         sRaw, [], 'use_closest_wl', 1);
                                     
            testCase.assertEmpty(sFilesOut);
            sInput = bst_process('GetInputStruct', sRaw);
            sStudy = bst_get('Study', sInput.iStudy);
            testCase.assertEqual(sStudy.iHeadModel, 1);
            testCase.assertEqual(sStudy.HeadModel.Comment, output_name);
            
            ilast_report = size(GlobalData.ProcessReports.Reports, 1);
            
             %% Call proc again
             sFilesOut = nst_run_bst_proc(output_name, 1, 'process_nst_import_head_model', ...
                                          sRaw, [], 'use_closest_wl', 1);
             testCase.assertEmpty(sFilesOut);
             testCase.assertTrue(~isempty(strfind(GlobalData.ProcessReports.Reports{ilast_report+1,4}, ...
                                                 'Force redo - removed previous result(s):')));
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

            resampled_data = in_bst_data(sFilesOut);
            testCase.assertMatches(resampled_data.Comment, output_name);
            testCase.assertTrue(all_close(resampled_data.Time, (0:2:(size(y,2)-1))*dt));
            
            %% Call proc again
            sFilesOut = nst_run_bst_proc('dummy_resampled', 0, 'process_resample', sRaw, [], 'freq', 5);
            
            testCase.assertTrue(exist(file_fullpath(sFilesOut), 'file')==2);
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end-1,2}.Comment, 'Set comment'); %from previous call
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,4}, ...
                                   'Skipped execution of process_resample. Outputs found.');                    
        end        
        
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

            testCase.assertTrue(isempty(sFilesOut));
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,1}, 'process');
            testCase.assertMatches(GlobalData.ProcessReports.Reports{end,2}.Comment, 'Delete files');
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
                                'Cannot safely manage unique outputs. Found duplicate items: dummy_raw/dummy_resampled')));

        end
        
        
        function test_run_proc_multiple_outputs(testCase)
            
        end
        
    end
end