classdef MontageTest < matlab.unittest.TestCase
     
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
 
        
        function test_unformat_channels(testCase)
            chan_types = nst_measure_types();
            
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
                testCase.assertMatches(ME.identifier, 'NIRSTORM:NonHomogeneousMeasure');
            end
            
            meas_nb_not_consistent = {'S01D03WL685', 'S01D03WL830', 'AUX1', 'S1D1WL685', 'S2D2WL830'};
            [isrcs, idets, measures, mtype] = nst_unformat_channels(meas_nb_not_consistent, 1);
            testCase.assertEqual(isrcs, [1 1 nan 1 2]);
            testCase.assertEqual(idets, [3 3 nan 1 2]);
            testCase.assertEqual(measures, [685 830 nan 685 830]);
            testCase.assertEqual(mtype, chan_types.WAVELENGTH);
            
            [warning_msg, warning_msg_id] = lastwarn;
            testCase.assertMatches(warning_msg, 'Inconsistent measure.*');
            
        end
  
        function test_montage_info(testCase)
            
            bst_create_test_subject();
            ref_point =  [-0.0789 ; -0.0263 ; 0.0655]; % meter, somewhere at the back of the head when using Colin27_4NIRS
            srcs_pos = [ref_point ref_point ref_point ref_point];
            dets_pos = [ref_point+[0.01;0;0] ref_point+[0.01;0;0] ref_point+[0.01;0.01;0.01] ref_point+[0.01;0.01;0.01]];
            nirs_input_file = bst_create_nirs_data('test', [], [], {'S2D2WL650', 'S2D2WL800', 'S2D3WL800', 'S2D3WL650'}, srcs_pos, dets_pos);
            
            sInput = in_bst_data(nirs_input_file);
            channel_def = in_bst_channel(sInput.ChannelFile);
            
            montage_info = nst_montage_info_from_bst_channels(channel_def.Channel); 
                       
            testCase.assertEqual(montage_info.src_ids, 2);
            testCase.assertTrue(all(sort(montage_info.det_ids) == [2 3]));
            
            testCase.assertTrue(all(montage_info.src_pos == ref_point'));
            testCase.assertTrue(all(all(montage_info.det_pos == [ref_point+[0.01;0;0] ref_point+[0.01;0.01;0.01]]')));
            
            testCase.assertTrue(all(montage_info.src_ichans{1} == [1 2 3 4]));
            
            testCase.assertTrue(all(sort(montage_info.det_ichans{find(montage_info.det_ids==2)}) == [1 2]));
            testCase.assertTrue(all(sort(montage_info.det_ichans{find(montage_info.det_ids==3)}) == [3 4]));
                
            for ipair=1:length(montage_info.pair_names)
               switch montage_info.pair_names{ipair}
                   case 'S2D2'
                       testCase.assertTrue(all(montage_info.pair_ichans(ipair,:) == [1 2])); 
                       testCase.assertTrue(all(montage_info.pair_sd_indexes(ipair,:) == [1 1]));
                       testCase.assertTrue(all(all(squeeze(montage_info.pair_loc(ipair,:,:)) == [srcs_pos(:,1) dets_pos(:,2)])));
                   case 'S2D3'
                       testCase.assertTrue(all(montage_info.pair_ichans(ipair,:) == [4 3])); 
                       testCase.assertTrue(all(montage_info.pair_sd_indexes(ipair,:) == [1 2]));
                       testCase.assertTrue(all(all(squeeze(montage_info.pair_loc(ipair,:,:)) == [srcs_pos(:,1) dets_pos(:,3)])));
                    otherwise
                       error(['Wrong pair:' montage_info.pair_names{ipair}]);
               end
            end
            
            %TODO: check and test process_nst_mbll/group_paired_channels
        end
        
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
            chan_types = nst_measure_types();
            
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
            
            bad_labels = {'S1D12Hb0', 'SS1D12HbR', 'S1D12wl685', 'AUX1', 'S1D12CCO'};
            for ilabel=1:length(bad_labels)
                [isrc, idet, measure, mtype] = nst_unformat_channel(bad_labels{ilabel}, 1);
                testCase.assertTrue(isnan(isrc));
                testCase.assertTrue(isnan(idet));
                testCase.assertTrue(isnan(measure));
                testCase.assertTrue(isnan(mtype));
                [warning_msg, warning_msg_id] = lastwarn;
                testCase.assertMatches(warning_msg, 'Malformed channel label.*');
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