classdef PrefixMatrixTest < matlab.unittest.TestCase
     
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

        function test_add_row_col_prefixes(testCase)
            
            ncols = 6;
            t1_col_names = arrayfun(@(n) sprintf('c%d', n), 1:ncols, 'UniformOutput', 0);
            t1_row_names = {'t1_r1', 't1_r2', 't1_r3'};
            t1 = array2table(randi(100, 3, ncols), 'VariableNames', t1_col_names, ...
                             'RowNames', t1_row_names);
                          
            subject_name = bst_create_test_subject('');
            
            sFile_t1 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table1');
            
            sFile_rt1 = bst_process('CallProcess', 'process_nst_prefix_matrix', ...
                                    sFile_t1, [], ...
                                    'row_prefixes', 'rp1_, rp2_, rp3_', 'col_prefixes', 'cp_');
            prefixed_mat = in_bst_matrix(sFile_rt1.FileName);
            for icol=1:length(prefixed_mat.ColNames)
                testCase.assertMatches(prefixed_mat.ColNames{icol}, ['cp_' t1_col_names{icol}]);
            end
            testCase.assertMatches(prefixed_mat.RowNames{1}, ['rp1_' t1_row_names{1}]);
            testCase.assertMatches(prefixed_mat.RowNames{2}, ['rp2_' t1_row_names{2}]);
            testCase.assertMatches(prefixed_mat.RowNames{3}, ['rp3_' t1_row_names{3}])
        end
        % TODO: test errors  
    end
end
