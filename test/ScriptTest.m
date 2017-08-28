classdef ScriptTest < matlab.unittest.TestCase
    
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
        function closeFigure(testCase)
            rmdir(testCase.tmp_dir, 's');
        end
    end
    
    methods(Test)
        function test_tutorial_tapping(testCase)
            addpath('../scripts');
            zip_fn = nst_request_files({{'tutorials', 'sample_nirs.zip'}}, ...
                                       1, 'ftp://neuroimage.usc.edu/pub/');
            assert(~isempty(zip_fn));
            unzip(zip_fn{1}, testCase.tmp_dir);
            nst_tutorial_tapping(testCase.tmp_dir, testCase.tmp_dir);
        end
    end
end

