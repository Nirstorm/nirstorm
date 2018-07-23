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
           
            delta_hb = [ 0.1;... %HbO (mmol.l-1)
                        -0.05];  %HbR (mmol.l-1)
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
%             figure();
%             plot(time,sDataOut.F);
            
            assert(all(all(abs(sDataOut.F(:,1:(activ_window_samples(1)-1))) < 1.37e-6))); %baseline goes close to 0
            
            % assert values are constant in activation window:
            assert(all(all(sDataOut.F(:,activ_window_samples) == sDataOut.F(:,activ_window_samples(1)))));
            
            % Check Hb values in activation window:
            hb_simu = hb_0 + delta_hb;
            baseline_mbll = sDataOut.F(1:2,1);
            hb_mbll = sDataOut.F(1:2,activ_window_samples(1)) - baseline_mbll;
            
            % non-regression test:
            non_reg_error = [0.6026e-4;0.9513e-4]; % quite high error...
            error = abs(hb_mbll - hb_simu);
            format short e;
            if all(error > non_reg_error)
                hb_table = table(hb_simu, hb_mbll, error);
                msg = sprintf(['Hb values from MBLL not consistent with simulated values:\n%s',...
                               evalc('disp(hb_table)')]);
                testCase.assertTrue(0, msg);
            elseif all( (error < non_reg_error) & (error > 1e-7) )
                hb_table = table(hb_simu, hb_mbll, error);
                msg = 'Hb values from MBLL not consistent with simulated values:\n%s';
                warning('Nirstorm:InaccurateResult', msg, evalc('disp(hb_table)'));
            end
            format short; % assume defatul short format was used
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

