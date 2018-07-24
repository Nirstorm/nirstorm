classdef SeparationTest < matlab.unittest.TestCase
    
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
        
        function test_forged_data(testCase)
            bst_create_test_subject();
            ref_point =  [-0.0789 ; -0.0263 ; 0.0655]; % meter, somewhere at the back of the head when using Colin27_4NIRS
            srcs_pos = [ref_point ref_point ref_point ref_point];
            dets_pos = [ref_point+[0.01;0;0] ref_point+[0.01;0;0] ref_point+[0.01;0.01;0.01] ref_point+[0.01;0.01;0.01]];
            nirs_input = bst_create_nirs_data('test', [], [], {}, srcs_pos, dets_pos);
            output = bst_process('CallProcess', 'process_nst_separations', nirs_input, []);
            
            sDataOut = in_bst_data(output.FileName);
            testCase.assertTrue(all(size(sDataOut.F)==[4 1]));
            testCase.assertTrue(all_close(sDataOut.F(1), 0.01 * 100)); %S1D1WL1
            testCase.assertTrue(all_close(sDataOut.F(2), 0.01 * 100)); %S1D1WL2
            expected_separation = sqrt(0.01^2 * 3) * 100;
            testCase.assertTrue(all_close(sDataOut.F(3), expected_separation)); %S1D2WL1
            testCase.assertTrue(all_close(sDataOut.F(4), expected_separation)); %S1D2WL2
        end
    end
end