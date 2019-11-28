function model = model_mc_splhblf_init_default()

model.var_hist_pace = 1;
model.obs_hist_pace = 1;
model.max_iterations = 100;
model.burnin_period = round(model.max_iterations/3);

model.options.response_duration = 25; % sec.
model.options.knots_type = 'denser_during_1st_half';
model.options.regular_knots_dt = 4; %sec.

model.options.trend_type = 'poly';
model.options.polybl_nb_coeffs = 4;

% For cosine baseline:
model.options.trend_type = 'cosine';
model.options.trend_bands(1).name = 'baseline';
model.options.trend_bands(1).bounds = [0 0.01];
model.options.trend_bands(2).name = 'respiration';
model.options.trend_bands(2).bounds = [0.16 0.33];
model.options.trend_bands(3).name = 'heart';
model.options.trend_bands(3).bounds = [1 1.25];
model.options.npb_mirror_trend_fix = 10;

% model.options.bandpass_residuals = [0.005 0.08]; %Hz

model.options.use_true_stim_induced = 0;
model.options.use_true_trend = 0;
model.options.use_true_noise = 0;

model.options.response_edge_tol = 1.5; %sec. tolerance to consider that ttp/ttu
                                       % is at the edge of the temporal
                                       % support -> used to clean the
                                       % posterior distribution and avoid
                                       % spurious bimodality if maximum is
                                       % located near the edge

var_names = {'f1', 'response', 'response_max', 'response_min', ...
             'response_time_to_peak', 'response_time_to_undershoot', ...
             'response_block', 'f1_var', 'noise_var', ...
             'trend_coeffs', 'trend_var', 'trend'};
         
for iv=1:length(var_names)
    model.variables_to_estimate.(var_names{iv}) = 1;
    model.variables.(var_names{iv}) = nan;
    model.check_rtol.(var_names{iv}) = 1e-5;
    model.check_atol.(var_names{iv}) = 1e-5;
end

end