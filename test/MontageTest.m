classdef MontageTest < matlab.unittest.TestCase
    methods(Test)
        
        function test_unformat_channel(testCase)
            chan_types = nst_channel_types();
            
            [isrc, idet, measure, mtype] = nst_unformat_channel('S01D03WL685');
            testCase.assertEqual(isrc, 1);
            testCase.assertEqual(idet, 3);
            testCase.assertEqual(measure, 685);
            testCase.assertEqual(mtype, chan_types.WAVELENGTH);
            
            [isrc, idet, measure, mtype] = nst_unformat_channel('S1D12HbT');
            testCase.assertEqual(isrc, 1);
            testCase.assertEqual(idet, 12);
            testCase.assertEqual(measure, 'HbT');
            testCase.assertEqual(mtype, chan_types.HB);
            
            bad_labels = {'S1D12Hb0', 'SS1D12HbR', 'S1D12wl685', 'S1D12CCO'};
            for ilabel=1:length(bad_labels)
                try
                    nst_unformat_channel(bad_labels{ilabel});
                catch ME
                    testCase.assertMatches(ME.identifier, 'NIRSTORM:MalformedChannelLabel');
                end
            end
        end

        function test_format_channel(testCase)
            testCase.assertTrue(strcmp(nst_format_channel(1, 3, 685), 'S1D3WL685'));
            
            testCase.assertTrue(strcmp(nst_format_channel(10, 14, 830), 'S10D14WL830'));
            
            testCase.assertTrue(strcmp(nst_format_channel(10, 14, 'HbO'),'S10D14HbO'));
            
            %TODO: test bad inputs
        end

        
    end
    
end