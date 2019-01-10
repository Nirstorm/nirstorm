classdef ProjectionTest < matlab.unittest.TestCase
    
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
        
        function test_two_channels(testCase)
            global GlobalData;
            [subject_name, sSubject] = bst_create_test_subject();
            
            dt = 0.05; %sec (20Hz)
            time = (0:8000) * dt; %sec
            nb_samples = length(time);
            
            
            % simulate activation signal
            delta_hb = [ 0.01;... %HbO (mmol.l-1)
                -0.005];  %HbR (mmol.l-1)
            delta_hb = delta_hb / 1000; % mol.l-1
            
            age = 25;
            % HbT_0 = 1.7 / 1000; % (mol.l-1) -> blood concentration
            HbT_0 = 100*10^-6; % (mol.l-1) -> tissue concentration
            sat = 0.7;
            hb_0 = [HbT_0 * sat;
                HbT_0 * (1-sat)];
            
            separation = 3; %cm
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
            
            y_dhb = zeros(2, nb_samples);
            y_dhb(:, activ_window_samples) = repmat(delta_hb, 1, length(activ_window_samples));
            
            signals = [y;y];
            nirs_input_mat_fn = bst_create_nirs_data('test_2_pairs', signals, time, ...
                {'S1D1WL690','S1D1WL832', ...
                'S1D2WL690','S1D2WL832'});
            
            % Generate fluences
            head_mesh_fn = sSubject.Surface(sSubject.iScalp).FileName;
            sMri = in_mri_bst(sSubject.Anatomy(sSubject.iAnatomy).FileName);
            sHead = in_tess_bst(head_mesh_fn);
            nirs_input = in_bst_data(nirs_input_mat_fn);
            ChannelMat = in_bst_channel(nirs_input.ChannelFile);
            montage_info = nst_montage_info_from_bst_channels(ChannelMat.Channel);
            src_coords = montage_info.src_pos;
            det_coords = montage_info.det_pos;
            [src_hv_idx det_hv_idx] = process_nst_import_head_model('get_head_vertices_closest_to_optodes', ...
                sMri, sHead, src_coords, det_coords);
            fluence_dir = cpt_spherical_fluences(sSubject, [src_hv_idx;det_hv_idx], ChannelMat.Nirs.Wavelengths, 50);
            
            for isrc=1:size(src_coords, 1)
                fprintf('S%d -> vertex #%d\n', isrc, src_hv_idx(isrc));
            end
            for idet=1:size(det_coords, 1)
                fprintf('D%d -> vertex #%d\n', idet, det_hv_idx(idet));
            end
            
            % Compute head model
            bst_process('CallProcess', ...
                'process_nst_import_head_model', nirs_input_mat_fn, [], ...
                'data_source', fluence_dir, ...
                'do_export_fluence_vol', 0, ...
                'outputdir', fluence_dir);
            
            % Compute dOD
            nirs_dOD = bst_process('CallProcess', ...
                'process_nst_dOD', nirs_input_mat_fn, [], ...
                'option_baseline_method', 1, ...  % mean
                'timewindow', [0 (activ_window_samples(1)-1)*dt]);
            
            % Do projection
            proj_methods = process_nst_cortical_projection('methods');
            result = bst_process('CallProcess', ...
                'process_nst_cortical_projection', nirs_dOD, [], ...
                'method', proj_methods.MNE);
            
            r1 = in_bst_results(result(1).FileName);
            head_model_1 = in_bst_headmodel(r1.HeadModelFile);
            [se_max, pos_max2] = max(squeeze(head_model_1.Gain(1,1,:)));
            r_hbo = r1.ImageGridAmp(pos_max2, :);
            
            r2 = in_bst_results(result(2).FileName);
            head_model_2 = in_bst_headmodel(r2.HeadModelFile);
            [se_max, pos_max2] = max(squeeze(head_model_1.Gain(1,1,:)));
            r_hbr = r2.ImageGridAmp(pos_max2, :);
            
            if 0
                % figure; hold on; plot(signals(1,:), 'k'); plot(signals(2,:), 'g');
                figure; hold on; plot(y_dhb(1,:), 'r'); plot(y_dhb(2,:), 'b');
                plot(r_hbo, 'r--');
                plot(r_hbr, 'b--');
                
            end
            
            % Non-regression tests
            % TODO: use tomographic forward model to generate input data and 
            % properly test results
            
            assert(abs(max(r_hbo) / max(abs(r_hbr)) - max(y_dhb(1,:)) / max(abs(y_dhb(2,:)))) < 1);
            
            assert(r_hbo(1) < 1e-6);
            assert(r_hbr(1) < 2e-7);
            
            assert(max(r_hbo) > 0);
            assert(min(r_hbr) < 0);

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
