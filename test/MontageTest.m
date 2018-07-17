classdef MontageTest < matlab.unittest.TestCase
    methods(Test)
        
        
        function test_pair_indexes(testCase)
            chan_names = {'S1D3WL685', 'S1D2WL685', 'S1D2WL830', 'S1D3WL830', ...
                          'S3D4WL830', 'S3D4WL685'};
            [pair_chan_indexes, pair_sd_ids] = nst_get_pair_indexes_from_names(chan_names);
            testCase.assertEqual(size(pair_chan_indexes, 1), 3);
            testCase.assertEqual(size(pair_chan_indexes, 2), 2);
            
            testCase.assertEqual(size(pair_sd_ids, 1), 3);
            testCase.assertEqual(size(pair_sd_ids, 2), 2);
            
            for ipair=1:size(pair_chan_indexes,1)
                switch(pair_chan_indexes(ipair, 1))
                    case 1
                       testCase.assertEqual(pair_chan_indexes(ipair, 2), 4);
                       testCase.assertEqual(pair_sd_ids(ipair,1), 1);
                       testCase.assertEqual(pair_sd_ids(ipair,2), 3);
                    case 2
                       testCase.assertEqual(pair_chan_indexes(ipair, 2), 3);
                       testCase.assertEqual(pair_sd_ids(ipair,1), 1);
                       testCase.assertEqual(pair_sd_ids(ipair,2), 2);
                    case 3
                       testCase.assertEqual(pair_chan_indexes(ipair, 2), 2);
                       testCase.assertEqual(pair_sd_ids(ipair,1), 1);
                       testCase.assertEqual(pair_sd_ids(ipair,2), 2);
                    case 4
                       testCase.assertEqual(pair_chan_indexes(ipair, 2), 1);
                       testCase.assertEqual(pair_sd_ids(ipair,1), 1);
                       testCase.assertEqual(pair_sd_ids(ipair,2), 3);
                    case 5
                       testCase.assertEqual(pair_chan_indexes(ipair, 2), 6);
                       testCase.assertEqual(pair_sd_ids(ipair,1), 3);
                       testCase.assertEqual(pair_sd_ids(ipair,2), 4);
                    case 6
                       testCase.assertEqual(pair_chan_indexes(ipair, 2), 5);
                       testCase.assertEqual(pair_sd_ids(ipair,1), 3);
                       testCase.assertEqual(pair_sd_ids(ipair,2), 4);
                    otherwise
                       testCase.assertTrue(0, sprintf('Wrong channel index %d', ...
                                                      pair_chan_indexes(ipair, 1))); 
                end
            end            
        end
        
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
                    testCase.assertTrue(0, 'NIRSTORM:MalformedChannelLabel not thrown');
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
        
        function test_unformat_channels(testCase)
            chan_types = nst_channel_types();
            
            good_chans_wl = {'S01D03WL685', 'S01D02WL685', 'S01D1WL685'};
            [isrcs, idets, measures, mtype] = nst_unformat_channels(good_chans_wl);
            testCase.assertEqual(isrcs, [1 1 1]);
            testCase.assertEqual(idets, [3 2 1]);
            testCase.assertEqual(measures, [685 685 685]);
            testCase.assertEqual(mtype, chan_types.WAVELENGTH);
            
            good_chans_hb = {'S01D02HbO', 'S1D2HbR', 'S1D1HbO', 'S1D1HbR'};
            [isrcs, idets, measures, mtype] = nst_unformat_channels(good_chans_hb);
            testCase.assertEqual(isrcs, [1 1 1 1]);
            testCase.assertEqual(idets, [2 2 1 1]);
            testCase.assertEqual(measures, {'HbO', 'HbR',  'HbO', 'HbR'});
            testCase.assertEqual(mtype, chan_types.HB);
            
            non_uniques = {'S01D03WL685', 'S1D3WL685', 'S01D1WL685'};
            try
                nst_unformat_channels(non_uniques);
                throw(MException('NIRSTORM:ExceptionNotThrown', 'Exception not thrown'));
            catch ME
                testCase.assertMatches(ME.identifier, 'NIRSTORM:NonUniqueChannels');
            end
            
            non_homogeneous = {'S01D03WL685', 'S1D2HbO', 'S01D1WL685'};
            try
                nst_unformat_channels(non_homogeneous);
                throw(MException('NIRSTORM:ExceptionNotThrown', 'Exception not thrown'));
            catch ME
                testCase.assertMatches(ME.identifier, 'NIRSTORM:NonHomogeneousChannelType');
            end
            
            meas_nb_not_consistent = {'S01D03WL685', 'S01D03WL830', 'S1D1WL685', 'S2D2WL830'};
            [isrcs, idets, measures, mtype] = nst_unformat_channels(meas_nb_not_consistent);
            testCase.assertEqual(isrcs, [1 1 1 2]);
            testCase.assertEqual(idets, [3 3 1 2]);
            testCase.assertEqual(measures, [685 830 685 830]);
            testCase.assertEqual(mtype, chan_types.WAVELENGTH);
            
            [warning_msg, warning_msg_id] = lastwarn;
            testCase.assertMatches(warning_msg, 'Inconsistent measure.*');
            
        end     
    end
end