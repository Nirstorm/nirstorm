classdef ConcatMatrixTest < matlab.unittest.TestCase
     
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

        function test_concatenate_by_row_errors(testCase)
            
            % Not the same columns
            ncols = 6;
            t1_col_names = arrayfun(@(n) sprintf('c%d', n), 1:ncols, 'UniformOutput', 0);
            t1_row_names = {'t1_r1', 't1_r2', 't1_r3'};
            t1 = array2table(randi(100, 3, ncols), 'VariableNames', t1_col_names, ...
                             'RowNames', t1_row_names);
                         
            t2_col_names = t1_col_names;
            t2_col_names{end} = 'waza';
            t2_row_names = {'t2_r1', 't2_r2', 't2_r3', 't2_r4', 't2_r5'};
            t2 = array2table(randi(100, 5, ncols), 'VariableNames', t2_col_names, ...
                             'RowNames', t2_row_names);
            
            subject_name = bst_create_test_subject('');
            
            sFile_t1 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table1');
            sFile_t2 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table2');
            
            stacking_types = process_nst_concat_matrices('get_stacking_types');
            try
                bst_process('CallProcess', 'process_nst_concat_matrices', ...
                            {sFile_t1, sFile_t2}, [], ...
                            'stacking_type', stacking_types.row);
                throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
            catch ME
                testCase.assertMatches(ME.identifier, 'Nirstorm:IncompatibleMatrices');
                testCase.assertMatches(ME.message, 'This file has columns incompatible with the first one.*');
            end            
        end
        
        function test_concatenate_by_col_errors(testCase)
            % Not the same rows
            nrows = 6;
            t1_row_names = arrayfun(@(n) sprintf('r%d', n), 1:nrows, 'UniformOutput', 0);
            t1_col_names = {'t1_c1', 't1_c2', 't1_c3'};
            t1 = array2table(randi(100, nrows, 3), 'VariableNames', t1_col_names, ...
                             'RowNames', t1_row_names);
            t2_col_names = {'t2_c1', 't2_c2', 't2_c3', 't2_c4', 't2_c5'};
            t2_row_names = t1_row_names;
            t2_row_names{end} = 'waza';
            t2 = array2table(randi(100, nrows, 5), 'VariableNames', t2_col_names, ...
                             'RowNames', t2_row_names);
            
            subject_name = bst_create_test_subject('');
            
            sFile_t1 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table1');
            sFile_t2 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table2');
            
            stacking_types = process_nst_concat_matrices('get_stacking_types');

            try
                bst_process('CallProcess', 'process_nst_concat_matrices', ...
                            {sFile_t1, sFile_t2}, [], ...
                            'stacking_type', stacking_types.column);
                throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
            catch ME
                testCase.assertMatches(ME.identifier, 'Nirstorm:IncompatibleMatrices');
                testCase.assertMatches(ME.message, 'This file has rows incompatible with the first one.*');
            end             
        end        
        
        function test_warnings(testCase)
            
            nrows = 6;
            row_names = arrayfun(@(n) sprintf('r%d', n), 1:nrows, 'UniformOutput', 0);
            t1_col_names = {'t1_c1', 't1_c2', 't1_c3'};
            t1 = array2table(randi(100, nrows, 3), 'VariableNames', t1_col_names, ...
                             'RowNames', row_names);
            t2_col_names = {'t2_c1', 't2_c2', 't2_c3', 't2_c4', 't2_c5'};
            t2 = array2table(randi(100, nrows, 5), 'VariableNames', t2_col_names, ...
                             'RowNames', row_names);
            
            subject_name = bst_create_test_subject('');    
                         
            extra.Events = db_template('Event');
            sFile_t1 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table1', extra);
            sFile_t2 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table2');
            
            stacking_types = process_nst_concat_matrices('get_stacking_types');
            bst_process('CallProcess', 'process_nst_concat_matrices', ...
                        {sFile_t1, sFile_t2}, [], ...
                        'stacking_type', stacking_types.column);
            
            
            [warning_msg, warning_msg_id] = lastwarn;
            testCase.assertMatches(warning_msg, 'Events field not empty.*');

            sFile_t3 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table3');
            sFile_t4 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table4', extra);
            
            bst_process('CallProcess', 'process_nst_concat_matrices', ...
                        {sFile_t3, sFile_t4}, [], ...
                        'stacking_type', stacking_types.column);
                     
            [warning_msg, warning_msg_id] = lastwarn;
            testCase.assertMatches(warning_msg, 'Events field not empty.*');
           
            
%             extra = struct();
%             extra.Time = 1:nrows;
%             sFile_t5 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table5', extra);
%             sFile_t6 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table6');
%             
%             bst_process('CallProcess', 'process_nst_concat_matrices', ...
%                         {sFile_t5, sFile_t6}, [], ...
%                         'stacking_type', stacking_types.column);
% 
%             [warning_msg, warning_msg_id] = lastwarn;
%             testCase.assertMatches(warning_msg, 'Time field not empty.*');
 
%             sFile_t7 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table7');
%             sFile_t8 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table8', extra);
%             
%             bst_process('CallProcess', 'process_nst_concat_matrices', ...
%                         {sFile_t7, sFile_t8}, [], ...
%                         'stacking_type', stacking_types.column);
%             
%             
%             [warning_msg, warning_msg_id] = lastwarn;
%             testCase.assertMatches(warning_msg, 'Time field not empty.*');
            
            
            
%             extra = struct();
%             extra.Description = t1_col_names;
%             sFile_t9 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table9', extra);
%             sFile_t10 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table10');
%             
%             bst_process('CallProcess', 'process_nst_concat_matrices', ...
%                         {sFile_t9, sFile_t10}, [], ...
%                         'stacking_type', stacking_types.column);
% 
%             [warning_msg, warning_msg_id] = lastwarn;
%             testCase.assertMatches(warning_msg, 'Description field not empty.*');
%  
%             extra.Description = t2_col_names;
%             sFile_t11 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table11');
%             sFile_t12 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table12', extra);
%             
%             bst_process('CallProcess', 'process_nst_concat_matrices', ...
%                         {sFile_t11, sFile_t12}, [], ...
%                         'stacking_type', stacking_types.column);
%             
%             
%             [warning_msg, warning_msg_id] = lastwarn;
%             testCase.assertMatches(warning_msg, 'Description field not empty.*');
            
        end
        
        
        function test_concatenate_by_row(testCase)
            ncols = 6;
            col_names = arrayfun(@(n) sprintf('c%d', n), 1:ncols, 'UniformOutput', 0);
            t1_row_names = {'t1_r1', 't1_r2', 't1_r3'};
            t1 = array2table(randi(100, 3, ncols), 'VariableNames', col_names, ...
                             'RowNames', t1_row_names);
            t2_row_names = {'t2_r1', 't2_r2', 't2_r3', 't2_r4', 't2_r5'};
            t2 = array2table(randi(100, 5, ncols), 'VariableNames', col_names, ...
                             'RowNames', t2_row_names);
            expected_merged_table = [t1 ; t2];
            expected_merged_values = table2array(expected_merged_table);
            expected_merged_row_names = [t1_row_names t2_row_names];
            
            subject_name = bst_create_test_subject('');
            
            sFile_t1 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table1');
            sFile_t2 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table2');
            
            stacking_types = process_nst_concat_matrices('get_stacking_types');

            sFiles_tmerged = bst_process('CallProcess', ...
                                         'process_nst_concat_matrices', ...
                                         {sFile_t1, sFile_t2}, [], ...
                                         'stacking_type', stacking_types.row);
            merged_mat = in_bst_matrix(sFiles_tmerged.FileName);
            
            testCase.assertTrue(all(merged_mat.Value(:) == expected_merged_values(:)));
            
            testCase.assertEqual(length(merged_mat.RowNames), length(expected_merged_row_names)); 
            for ir=1:length(expected_merged_row_names)
                testCase.assertMatches(expected_merged_row_names{ir}, merged_mat.RowNames{ir});
            end
            
            testCase.assertEqual(length(merged_mat.ColNames), length(col_names));
            for ic=1:length(col_names)
                testCase.assertMatches(col_names{ic}, merged_mat.ColNames{ic});
            end
        end
        
        function test_concatenate_by_col(testCase)
            nrows = 6;
            row_names = arrayfun(@(n) sprintf('r%d', n), 1:nrows, 'UniformOutput', 0);
            t1_col_names = {'t1_c1', 't1_c2', 't1_c3'};
            t1 = array2table(randi(100, nrows, 3), 'VariableNames', t1_col_names, ...
                             'RowNames', row_names);
            t2_col_names = {'t2_c1', 't2_c2', 't2_c3', 't2_c4', 't2_c5'};
            t2 = array2table(randi(100, nrows, 5), 'VariableNames', t2_col_names, ...
                             'RowNames', row_names);
            expected_merged_table = [t1 t2];
            expected_merged_values = table2array(expected_merged_table);
            expected_merged_col_names = [t1_col_names t2_col_names];
            
            subject_name = bst_create_test_subject('');
            
            sFile_t1 = nst_save_table_in_bst(t1, subject_name, 'condition', 'table1');
            sFile_t2 = nst_save_table_in_bst(t2, subject_name, 'condition', 'table2');
            
            stacking_types = process_nst_concat_matrices('get_stacking_types');

            sFiles_tmerged = bst_process('CallProcess', ...
                                         'process_nst_concat_matrices', ...
                                         {sFile_t1, sFile_t2}, [], ...
                                         'stacking_type', stacking_types.column);
            merged_mat = in_bst_matrix(sFiles_tmerged.FileName);
            
            testCase.assertTrue(all(merged_mat.Value(:) == expected_merged_values(:)));
            
            testCase.assertEqual(length(merged_mat.RowNames), length(row_names)); 
            for ir=1:length(row_names)
                testCase.assertMatches(row_names{ir}, merged_mat.RowNames{ir});
            end
            
            testCase.assertEqual(length(merged_mat.ColNames), length(expected_merged_col_names));
            for ic=1:length(expected_merged_col_names)
                testCase.assertMatches(expected_merged_col_names{ic}, merged_mat.ColNames{ic});
            end
        end        
    end
end
