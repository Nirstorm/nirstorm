classdef FilterTableTest < matlab.unittest.TestCase
    
    properties
        acq_table
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            s_table(1).subject_id = 1;
            s_table(1).period = 'pre';
            s_table(1).acq_date = '2017/05/23';
            s_table(1).age = 65;
            s_table(1).task = 'NB';
            s_table(1).data_fn = '/path/to/2017_05_23_NB.nirs';
            
            s_table(2).subject_id = 1;
            s_table(2).period = 'pre';
            s_table(2).acq_date = '2017/05/23';
            s_table(2).age = 65;
            s_table(2).task = 'DT';
            s_table(2).data_fn = '/path/to/2017_05_23_DT.nirs';
            
            s_table(3).subject_id = 1;
            s_table(3).period = 'mid';
            s_table(3).acq_date = '2017/08/22';
            s_table(3).age = 65;
            s_table(3).task = 'NB';
            s_table(3).data_fn = '/path/to/2017_08_22_NB.nirs';
            
            s_table(4).subject_id = 2;
            s_table(4).period = 'pre';
            s_table(4).acq_date = '2016/12/01';
            s_table(4).age = 42;
            s_table(4).task = 'NB';
            s_table(4).data_fn = '/path/to/2016_12_01_NB.nirs';
            
            s_table(5).subject_id = 2;
            s_table(5).period = 'pre';
            s_table(5).acq_date = '2016/12/01';
            s_table(5).age = 42;
            s_table(5).task = 'DT';
            s_table(5).data_fn = '/path/to/2016_12_01_DT.nirs';
         
            s_table(6).subject_id = 2;
            s_table(6).period = 'post';
            s_table(6).acq_date = '2017/09/22';
            s_table(6).age = 43;
            s_table(6).task = 'NB';
            s_table(6).data_fn = '/path/to/2017_09_22_NB.nirs';
            
            testCase.acq_table = struct2table(s_table);
        end
    end    
    
        methods(Test)
            
            function test_include_values(testCase)
                filters_acquisitions.period.include = {'pre', 'post'};
                filtered_table = nst_filter_table(testCase.acq_table, filters_acquisitions);
                
                testCase.assertEqual(size(filtered_table, 1), 5);
                testCase.assertTrue(all(ismember(filtered_table.period), {'pre', 'post'}));
            end
            
            function test_include_function(testCase)
                filters_acquisitions.age.include = @(age) age < 60 ;
                filtered_table = nst_filter_table(testCase.acq_table, filters_acquisitions);
                
                testCase.assertEqual(size(filtered_table, 1), 3);
                testCase.assertEqual(unique(filtered_table.subject_id), 1);  
            end
            
            function test_exclude_values(testCase)
                filters_acquisitions.subject_id.exclude = 1;
                filtered_table = nst_filter_table(testCase.acq_table, filters_acquisitions);
                
                testCase.assertEqual(size(filtered_table, 1), 3);
                testCase.assertEqual(unique(filtered_table.subject_id), 2);    
            end
            
            function test_exclude_function(testCase)
                filters_acquisitions.age.include = @(age) age < 60 ;
                filtered_table = nst_filter_table(testCase.acq_table, filters_acquisitions);
                
                testCase.assertEqual(size(filtered_table, 1), 3);
                testCase.assertEqual(unique(filtered_table.subject_id), 2);  
            end
            
            function test_entry_include(testCase)
                filters_acquisitions.entry__.include = [...
                          struct('subject_id', 1, 'period', 'pre' , 'task', 'NB'), ...
                          struct('subject_id', 2, 'period', 'pre' , 'task', 'DT'), ...
                          struct('subject_id', 10, 'period', 'dummy' , 'task', 'dummy')];
                
                filtered_table = nst_filter_table(testCase.acq_table, filters_acquisitions);
                
                testCase.assertEqual(size(filtered_table, 1), 2);
                testCase.assertEqual(filtered_table.subject_id(1), 1);
                testCase.assertEqual(filtered_table.subject_id(2), 2);
                testCase.assertEqual(filtered_table.period{1}, 'pre');
                testCase.assertEqual(filtered_table.period{2}, 'mid');
                testCase.assertEqual(filtered_table.task{1}, 'NB');
                testCase.assertEqual(filtered_table.task{2}, 'DT');
                testCase.assertEqual(filtered_table.acq_date{1}, '2017/05/23');
                testCase.assertEqual(filtered_table.acq_date{2}, '2016/12/01');
            end
            
            function test_entry_exclude(testCase)
                filters_acquisitions.entry__.exclude = [...
                          struct('subject_id', 1, 'period', 'pre' , 'task', 'NB'), ...
                          struct('subject_id', 2, 'period', 'pre' , 'task', 'DT'), ...
                          struct('subject_id', 10, 'period', 'dummy' , 'task', 'dummy')];
                
                filtered_table = nst_filter_table(testCase.acq_table, filters_acquisitions);
                
                testCase.assertEqual(size(filtered_table, 1), 4);
                testCase.assertEqual(filtered_table.subject_id(1), 1);
                testCase.assertEqual(filtered_table.subject_id(2), 1);
                testCase.assertEqual(filtered_table.subject_id(end), 2);
                
                testCase.assertEqual(filtered_table.period{1}, 'pre');
                testCase.assertEqual(filtered_table.period{2}, 'mid');
                testCase.assertEqual(filtered_table.period{end}, 'post');

                testCase.assertEqual(filtered_table.task{1}, 'DT');
                testCase.assertEqual(filtered_table.task{2}, 'NB');
                testCase.assertEqual(filtered_table.task{end}, 'NB');
                
                testCase.assertEqual(filtered_table.acq_date{1}, '2017/05/23');
                testCase.assertEqual(filtered_table.acq_date{2}, '2017/08/22');
                testCase.assertEqual(filtered_table.acq_date{end}, '2017/09/22');
            end
            
            
            function test_filter_invalid_column(testCase)
                filters_acquisitions.bad_col.include = 1;
                try
                    nst_filter_table(testCase.acq_table, filters_acquisitions);
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:InvalidFilterColumn');
                end
            
            end
            
            function test_filter_invalid_operation(testCase)
                filters_acquisitions.age.bad_operation = 1;
                try
                    nst_filter_table(testCase.acq_table, filters_acquisitions);
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:InvalidFilterOperation');
                end
            
            end
            
            function test_filter_invalid_value(testCase)
                filters_acquisitions.age.include = struct('yo', 1);
                try
                    nst_filter_table(testCase.acq_table, filters_acquisitions);
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:InvalidFilterOperation');
                end
             end
            
            function test_filter_inconsistent_value_type(testCase)
                filters_acquisitions.age.include = '1';
                try
                    nst_filter_table(testCase.acq_table, filters_acquisitions);
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:InconsistentFilterValue');
                end
            end
             
            function test_filter_function_not_predicate(testCase)
                filters_acquisitions.age.include = @(age) 'too_old';
                try
                    nst_filter_table(testCase.acq_table, filters_acquisitions);
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:FilterFunctionNotPredicate');
                end
            end
            
            function test_entry__in_talbe(testCase)
                entry__ = {'a', 'b', 'c'}';
                subject_id = (1:3)';
                atable = table(entry__, subject_id);
                try
                    nst_filter_table(atable, struct());
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:UnsupportedTableWithEntry__');
                end
            end

            function test_non_unique_entry(testCase)
                filters_acquisitions.entry__.include =[struct('subject_id', 1, 'period', 'pre')];
                try
                    nst_filter_table(testCase.acq_table, filters_acquisitions);
                    throw(MException('Nirstorm:ExceptionNotThrown', 'Exception not thrown'));
                catch ME
                    testCase.assertMatches(ME.identifier, 'Nirstorm:NonUniqueEntry');
                end
            end
            
        end
end