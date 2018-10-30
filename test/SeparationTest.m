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
        
        
        function test_forged_data_by_pairs(testCase)
            bst_create_test_subject();
            ref_point =  [-0.0789 ; -0.0263 ; 0.0655]; % meter, somewhere at the back of the head when using Colin27_4NIRS
            srcs_pos = [ref_point ref_point ref_point ref_point];
            dets_pos = [ref_point+[0.01;0;0] ref_point+[0.01;0;0] ref_point+[0.01;0.01;0.01] ref_point+[0.01;0.01;0.01]];
            nirs_input_file = bst_create_nirs_data('test', [], [], {'S2D2WL650', 'S2D2WL800', 'S2D1WL650', 'S2D1WL800'}, srcs_pos, dets_pos);
            
            pair_ids = [[2 1];[2 2]];
            sInput = in_bst_data(nirs_input_file);
            channel_def = in_bst_channel(sInput.ChannelFile);
            separations = process_nst_separations('Compute', channel_def.Channel, pair_ids);
            
            testCase.assertTrue(all(size(separations)==[size(pair_ids, 1) 1]));
            expected_separation = sqrt(0.01^2 * 3) * 100;
            testCase.assertTrue(all_close(separations(1), expected_separation)); %S2D1
            expected_separation = 0.01 * 100;
            testCase.assertTrue(all_close(separations(2), expected_separation)); %S2D2
        end
        
        function test_forged_data_by_channels(testCase)
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