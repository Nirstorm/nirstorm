classdef MbllTest < matlab.unittest.TestCase
    
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
        
        function test_mbll_with_simulations(testCase)
           
            delta_hb = [ 0.01;... %HbO (mmol.l-1)
                        -0.005];  %HbR (mmol.l-1)
            delta_hb = delta_hb / 1000; % mol.l-1
            
            age = 25;
            % HbT_0 = 1.7 / 1000; % (mol.l-1) -> blood concentration
            HbT_0 = 100*10^-6; % (mol.l-1) -> tissue concentration
            sat = 0.7;
            hb_0 = [HbT_0 * sat;
                    HbT_0 * (1-sat)];

            separation = 1; %cm
            wavelengths = [690, 832];
            pvf = 1;
            hb_extinctions = process_nst_mbll('get_hb_extinctions', wavelengths);
            i_light_ref = 1e6;
            
            y_baseline = mbll_fwd(i_light_ref, hb_0, hb_extinctions, age, separation, pvf);
            y_activ = mbll_fwd(i_light_ref, hb_0 + delta_hb, hb_extinctions, age, separation, pvf);
            
            dt = 0.1; %sec
            nb_samples = 1000;
            time = (0:(nb_samples-1))*dt;
            
            y = zeros(2, nb_samples) + repmat(y_baseline, 1, nb_samples);
            activ_window_samples = 300:500;
            y(:, activ_window_samples) = repmat(y_activ, 1, length(activ_window_samples));
            
            bst_create_test_subject();
            sRaw = bst_create_nirs_data('simulated_evoked_raw', y, time, {'S1D1WL690','S1D1WL832'},...
                                        [1 1;1 1;1 1], [1.01 1.01;1 1;1 1]);

            sOutput = bst_process('CallProcess', 'process_nst_mbll', sRaw, [], ...
                                  'option_age',             age, ...
                                  'option_pvf',             pvf, ...
                                  'option_baseline_method', 1, ...  % mean
                                  'timewindow', [0 (activ_window_samples(1)-1)*dt], ...
                                  'option_do_plp_corr',     1);
            sDataOut = in_bst_data(sOutput.FileName);
            
            dhb_signal = zeros(3, nb_samples);
            dhb_signal(1:2, activ_window_samples) = repmat(delta_hb, 1, length(activ_window_samples));
            dhb_signal(3,:) = sum(dhb_signal(1:2, :));
            
            if 0
                figure(); hold on; 
                plot(time,dhb_signal(1,:), 'r');
                plot(time,dhb_signal(2,:), 'b');
                
                plot(time, sDataOut.F(1,:), 'r--');
                plot(time, sDataOut.F(2,:), 'b--');
            end
            
            assert(all_close(dhb_signal, sDataOut.F));            
        end
        
        function test_mbll_on_tapping_data(testCase)
            repo_url = nst_get_repository_url();
            data_fns = nst_request_files({{'unittest','motor_data','motor.nirs'}, ...
                                          {'unittest','motor_data','optodes.txt'}, ...
                                          {'unittest','motor_data','fiducials.txt'}...
                                          {'unittest','motor_data','mbll.mat'}}, ...
                                         1, repo_url);
                                    
            nirs_fn = data_fns{1};
            sFile = utest_import_nirs_in_bst(nirs_fn);
            sFile_bad_tagged = bst_process('CallProcess', 'process_nst_detect_bad', sFile, [], ...
                                       'option_remove_negative', 1, ...
                                       'option_invalidate_paired_channels', 1, ...
                                       'option_max_sat_prop', 1);
            sFileHb = bst_process('CallProcess', 'process_nst_mbll', sFile_bad_tagged, [], ...
                        'option_age', {25, ''}, ...
                        'option_pvf', {50, ''}, ...
                        'option_baseline_method', {1, {'mean', 'median'}}, ...
                        'option_do_plp_corr', 1);
                    
            sDataOut = in_bst_data(sFileHb.FileName);
            
            % Non-regression test (on 23rd July 2018)
            expected_mbll_fn = data_fns{end};
            expected_mbll = load(expected_mbll_fn);
            expected_mbll = expected_mbll.mbll_motor;
            
            if 0
                figure(); hold on; 
                plot(sDataOut.Time,expected_mbll(1,:), 'r');
                plot(sDataOut.Time,expected_mbll(2,:), 'b');
                
                plot(sDataOut.Time, sDataOut.F(1,:), 'r--');
                plot(sDataOut.Time, sDataOut.F(2,:), 'b--');
            end
            
            testCase.assertTrue(all_close(sDataOut.F, expected_mbll));
        end
        
    end
    
end

function i_light_output = mbll_fwd(i_light_ref, concentrations, extinctions, age, separation, pvf)
y0 = [5.38;4.67];
a1 = [0.049;0.062];
a2 = [0.877;0.819];
dpf = y0 + a1 .* age.^a2; %checked that dfp was OK
i_light_output = i_light_ref * power(repmat(10,2,1), -separation .* extinctions * concentrations .* dpf / pvf);
end

