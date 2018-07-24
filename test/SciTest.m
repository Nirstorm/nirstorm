classdef SciTest < matlab.unittest.TestCase
    
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
        
        function test_sci_with_simulations(testCase)
            
            dt = 0.05; %sec (20Hz)
            time = (0:8000) * dt; %sec
            nb_samples = length(time);
            
            bst_create_test_subject();
            
            %% Only noise -- uncorrelated wavelength signals
            signals = [randn(1, nb_samples);randn(1, nb_samples);];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI close to 0
            testCase.assertTrue(all(sDataOut.F<0.1));
            
            %% Correlation outside of cardiac band
            freq = 5; % Hz
            signals = [cos(4*pi*freq*time) + randn(1, nb_samples);...
                       cos(4*pi*freq*time) + randn(1, nb_samples)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI close to 0
            testCase.assertTrue(all(sDataOut.F<0.1));
            
            %% Perfect cardiac consistency - correlated
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time); cos(4*pi*freq*time)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI high
            testCase.assertTrue(all(sDataOut.F>0.99));
            
            %% Perfect cardiac consistency - anticorrelated
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time); -cos(4*pi*freq*time)];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI high
            testCase.assertTrue(all(sDataOut.F>0.99));
                        
            %% Consistent cardiac components - high SNR
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time) + randn(1, nb_samples)*0.2; ...
                       -cos(4*pi*freq*time) + randn(1, nb_samples)*0.2];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI > 0.75
            testCase.assertTrue(all(sDataOut.F>0.75));
            
            %% Consistent cardiac components - very low SNR
            freq = 1; % Hz
            signals = [cos(4*pi*freq*time) + randn(1, nb_samples)*50; ...
                       -cos(4*pi*freq*time) + randn(1, nb_samples)*50];
            nirs_input = bst_create_nirs_data('test_all_noisy', signals, time);
            output = bst_process('CallProcess', 'process_nst_sci', nirs_input, []);
            sDataOut = in_bst_data(output.FileName);
            
            %Check SCI < 0.75
            testCase.assertTrue(all(sDataOut.F<0.75));
            
            
        end
        
        function test_sci_on_tapping_data(testCase)
            repo_url = nst_get_repository_url();
            data_fns = nst_request_files({{'unittest','motor_data','motor.nirs'}, ...
                                          {'unittest','motor_data','optodes.txt'}, ...
                                          {'unittest','motor_data','fiducials.txt'},...
                                          {'unittest','motor_data','sci.mat'}}, ...
                                         1, repo_url);

                                    
            nirs_fn = data_fns{1};
            sFile = utest_import_nirs_in_bst(nirs_fn);
            output = bst_process('CallProcess', 'process_nst_sci', sFile, []);
            sDataOut = in_bst_data(output.FileName);
            
            % Non-regression test (sci results manually checked and
            % recorded on 17 July 2018)
            expected_sci_fn = data_fns{end};
            expected_sci = load(expected_sci_fn);
            expected_sci = expected_sci.sci_motor;
            
            testCase.assertTrue(all_close(sDataOut.F, expected_sci, 0.01, 0.01));
        end
        
        
    end
end
