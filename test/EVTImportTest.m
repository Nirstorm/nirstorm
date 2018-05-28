classdef EVTImportTest < matlab.unittest.TestCase
    
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
        function test_csv_event_import_time_unit(testCase)
        end
    end
end

function events = get_events(sFile)
DataMat = in_bst_data(sFile.FileName, 'F');
events = DataMat.F.events;
end
