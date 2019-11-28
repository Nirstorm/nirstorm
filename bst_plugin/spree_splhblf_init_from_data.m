function  [model] = spree_splhblf_init_from_data(model, nirs)
% Args:
%   - data: NirsData (see nirs_data.m)
%

model = init_input_signals(model, nirs);
model = init_variables(model, nirs);

end

function model = init_input_signals(model, nirs)
% 
% if length(nirs.measure_labels) > 1
%     throw(MException('NirsDyna:WrongInputData', ...
%                      ['Input Nirs data should only contain one measure '...
%                       'e.g. only 690nm or only HbO']))
% end

model.constants.y = nirs.data;
model.constants.y_time_axis = nirs.time;
model.constants.dt = diff(nirs.time(1:2));
model.constants.n_channels = size(nirs.data, 2);
model.constants.n_samples = size(nirs.data, 1);
model.constants.sampling_freq = 1/model.constants.dt;
model.constants.npb_mirror_trend_fix = model.options.npb_mirror_trend_fix;
end

function model = init_variables(model, nirs)

if isfield(nirs, 'simulation')
    simulation = nirs.simulation;
else
    simulation = struct();
end

dt = 1/model.constants.sampling_freq;

% Init response basis matrix
if isfield(simulation, 'F')
    model.constants.F = simulation.F;
    model.constants.response_time_axis = simulation.response_time_axis;
    response_duration = length(model.constants.response_time_axis);
else
    
    if isfield(simulation, 'response_time_axis')
        model.constants.response_time_axis = simulation.response_time_axis;
        response_duration = simulation.response_time_axis(end); %sec
    else
        response_duration = model.options.response_duration; %sec
        model.constants.response_time_axis = (0:dt:response_duration)';
    end
    
    spline_order = 3;
    if strcmp(model.options.knots_type, 'regular')
        knots_dt = model.options.regular_knots_dt;
        knots = 0:knots_dt:response_duration;
    elseif strcmp(model.options.knots_type, 'denser_during_1st_half')
        assert(response_duration > 20);
        knots = [0:3:12 19 response_duration];
        %knots = [0:4:response_duration];
    else
        throw(MException('NIRSDyna:ScenarioNotImplemented'));
    end
    %knots = augknt(knots, floor(spline_order/2) + 1);
    
    if 1
        % Add double knots at start and end
        knots = [knots(1) knots response_duration];
    end

    model.constants.F = spcol(knots, spline_order, ...
        model.constants.response_time_axis);
    if 0
        fig_id = 2;
        fig_visibility = 'on';
        h = figure(fig_id); clf(fig_id); set(h, 'visible', fig_visibility); hold on;
        plot(model.constants.response_time_axis, model.constants.F);
        xlim(model.constants.response_time_axis([1 end]));
    end
end
model.constants.phi_nb_coeffs = size(model.constants.F, 2);
model.constants.response_nb_coeffs = length(model.constants.response_time_axis);

if isfield(nirs, 'simulation') && isfield(nirs.simulation, 'X')
    model.constants.X = nirs.simulation.X;
    assert(size(model.constants.X, 1) == model.constants.n_samples)
elseif isfield(nirs, 'paradigm')
    session_duration = model.constants.y_time_axis(end);
    if length(nirs.paradigm.onsets) > 1
        warning('Only 1st experimental conditions is used (TODO: handle multi-cond)');
    end
    model.constants.X = nst_build_onset_toeplitz_mtx(nirs.paradigm.onsets{1}, ...
                                                     nirs.paradigm.durations{1}, ...
                                                     1/model.constants.sampling_freq, ...
                                                     session_duration, ...
                                                     model.constants.response_nb_coeffs);
    assert(size(model.constants.X, 1) >= size(model.constants.y, 1));
    model.constants.X = model.constants.X(1:size(model.constants.y, 1), :);
    model.constants.paradigm = nirs.paradigm;
end
block_duration = nirs.paradigm.durations{1}(1) + model.constants.response_time_axis(end);
model.constants.X_block = nst_build_onset_toeplitz_mtx([0], ...
    nirs.paradigm.durations{1}(1), ...
    1/model.constants.sampling_freq, ...
    block_duration, ...
    model.constants.response_nb_coeffs);
model.constants.response_block_time_axis = (1:size(model.constants.X_block,1))' * dt;

% adjust threshold for proba
% 1% of variation within expected hemodynamics band
%TODO: expose options
response_length = model.constants.response_nb_coeffs;
response_cano = spm_hrf(dt, [5., 16, 1, 1, 6, 0, response_duration]);
if length(response_cano) > response_length
    response_cano = response_cano(1:response_length);
elseif length(response_cano) < response_length
    response_cano = [response_cano ; zeros(response_length-length(response_cano))];
end
    
expected_hemodyn = model.constants.X * response_cano;
% fig_from_label('spectrum_expected_hemodyn');
[pxx_hemodyn, f] = periodogram(expected_hemodyn,rectwin(length(expected_hemodyn)), ...
                               length(expected_hemodyn),model.constants.sampling_freq);
[spxx, idx_f_sorted] = sort(pxx_hemodyn, 'descend');

if 0
    fig_from_label('spectrum_expected_hemodyn');
    periodogram(expected_hemodyn,rectwin(length(expected_hemodyn)), ...
        length(expected_hemodyn),model.constants.sampling_freq);
    fig_from_label('expected_hemodyn');
    plot(model.constants.y_time_axis, expected_hemodyn);
    fig_from_label('filtered expected_hemodyn');
    opts.ftype = 'bandpass';
    if idx_f_sorted(2) ~= 1
        sf1 = f(idx_f_sorted(2));
    else
        sf1 = f(idx_f_sorted(1));
    end
    sf2 = f(idx_f_sorted(3));
    opts.low_cutoff = min(sf1, sf2);
    opts.high_cutoff = max(sf1, sf2);
    opts.order = 3;
    opts.flag_keepMean = 0;
    filt_expected_hemodyn = mfip_IIR_butter(expected_hemodyn, model.constants.sampling_freq , opts);
    plot(model.constants.y_time_axis, filt_expected_hemodyn);
end

opts.ftype = 'bandpass';
if ~(isfield(model.options, 'trend_bands') && any(strcmp('baseline', {model.options.trend_bands.name})))
    if idx_f_sorted(2) ~= 1
        sf1 = f(idx_f_sorted(2));
    else
        sf1 = f(idx_f_sorted(1));
    end
    sf2 = f(idx_f_sorted(3));
else
    sf1 = model.options.trend_bands(strcmp('baseline', {model.options.trend_bands.name})).bounds(2);
    sf2 = 0.3; %TODO: maybe use less arbitrary value
end
opts.low_cutoff = min(sf1, sf2);
opts.high_cutoff = max(sf1, sf2);
opts.order = 3;
opts.flag_keepMean = 0;
y_hemodyn = mfip_IIR_butter(model.constants.y, model.constants.sampling_freq , opts);
y_hemodyn = detrend(y_hemodyn);

model.summary.pk_threshold = max(y_hemodyn) * 0.0075; %0.01 %TODO: expose option

% model.constants.XtY = model.constants.X' * model.constants.y;
% model.constants.XtX = model.constants.X' * model.constants.X;
% model.constants.FFt = model.constants.F * model.constants.F';
model.constants.FtXt = model.constants.F' * model.constants.X';
model.constants.FtXtXF = model.constants.F' * (model.constants.X' * model.constants.X) * model.constants.F;
model.constants.XF = model.constants.X * model.constants.F;
%model.constants.FFt = model.constants.F * model.constants.F';

if isfield(simulation, 'P1')
    P1 = simulation.P1;
    D = simulation.D;
else
    [P1, D] = eigs(model.constants.F * model.constants.F', ...
                   size(model.constants.F, 2));
end
model.constants.D = D;
model.constants.P1 = P1;
model.constants.f1_nb_coeffs = size(P1, 2);
model.constants.XP = model.constants.X * P1;
model.constants.XtX = model.constants.X' * model.constants.X;
model.constants.PtXt = P1' * model.constants.X';
model.constants.PXXP = P1' * model.constants.XtX * P1;
% model.constants.PtXtY = P1' * model.constants.XtY; 

% f1 (response coeff)
response_init_base = [0 ; 0.05; 0.05; zeros(model.constants.response_nb_coeffs-3, 1)];
f1_init_base = model.constants.P1' * response_init_base;
f1_init = repmat(f1_init_base, 1, model.constants.n_channels);
if isfield(simulation, 'f1')
    model = init_variable(model, 'f1', {'response_coeff', 'channel'}, ...
                          f1_init, simulation.f1);
else
    model = init_variable(model, 'f1', {'response_coeff', 'channel'}, ...
                          f1_init);
end

% response
response_init = model.constants.P1 * f1_init;
if isfield(simulation, 'response')
    model = init_variable(model, 'response', {'response_time', 'channel'}, ...
                          response_init, simulation.response);
else
    model = init_variable(model, 'response', {'response_time', 'channel'}, ...
                          response_init);
end

response_block_init = model.constants.X_block * response_init;
model = init_variable(model, 'response_block', {'response_time', 'channel'}, ...
                      response_block_init);

% response peak
[init_peak, ittp] = max(response_init', [], 2);
init_time_to_peak = model.constants.response_time_axis(ittp);
if isfield(simulation, 'response')
    [sim_peak_value, ittp] = max(simulation.response);
    sim_ttp = simulation.response_time_axis(ittp);
    model = init_variable(model, 'response_max', {'channel'}, ...
                          init_peak, sim_peak_value');
    model = init_variable(model, 'response_time_to_peak', {'channel'}, ...
                          init_time_to_peak, sim_ttp);
else
    model = init_variable(model, 'response_max', {'channel'}, ...
                          init_peak');
    model = init_variable(model, 'response_time_to_peak', {'channel'}, ...
                          init_time_to_peak);
end

% response undershoot
[init_undershoot, ittu] = min(response_init', [], 2);
init_time_to_undershoot = model.constants.response_time_axis(ittu);
if isfield(simulation, 'response')
    [sim_undershoot_value, ittu] = min(simulation.response);
    sim_ttu = simulation.response_time_axis(ittu);
    model = init_variable(model, 'response_min', {'channel'}, ...
                          init_undershoot, sim_undershoot_value');
    model = init_variable(model, 'response_time_to_undershoot', {'channel'}, ...
                          init_time_to_undershoot, sim_ttu);
else
    model = init_variable(model, 'response_min', {'channel'}, ...
                          init_undershoot);
    model = init_variable(model, 'response_time_to_undershoot', {'channel'}, ...
                          init_time_to_undershoot);
end


% f1 variance
if isfield(simulation, 'f1_var')
    model = init_variable(model, 'f1_var', {'channel'}, ...
                          var(f1_init), simulation.f1_var);
else
    model = init_variable(model, 'f1_var', {'channel'}, ...
                          var(f1_init));
end

% noise variance
if isfield(simulation, 'noise_var_emp')
    model = init_variable(model, 'noise_var', {'channel'}, ...
                          ones(1, model.constants.n_channels) * .1,  ...
                          simulation.noise_var_emp);
else
    model = init_variable(model, 'noise_var', {'channel'}, ...
                          var(model.constants.y, 1) * .5);
end


% trend coeffs
if isfield(simulation, 'trend_coeffs')
    trend_bsize = size(simulation.trend_coeffs, 1);
    model = init_variable(model, 'trend_coeffs', {'coeff', 'channel'}, ...
                          randn(trend_bsize, model.constants.n_channels) * .1,  ...
                          simulation.trend_coeffs);
    model.constants.T = simulation.trend_mat;
    if isfield(simulation, 'npb_mirror_trend_fix')
        model.constants.npb_mirror_trend_fix = simulation.npb_mirror_trend_fix;
    end
else
    n_samples_trend = model.constants.n_samples + 2 * model.constants.npb_mirror_trend_fix;
    if strcmp(model.options.trend_type, 'poly')
        trend_bsize = model.options.polybl_nb_coeffs;
        trend_mat = build_basis_poly(n_samples_trend, trend_bsize, dt);
        model.constants.trend_bands.name = 'baseline';
        model.constants.trend_bands.bounds = [];
        model.constants.trend_bands.coeff_indexes = 1:trend_bsize;
    elseif strcmp(model.options.trend_type, 'cosine')
        [trend_mat, band_indexes] = build_basis_dct(n_samples_trend, 1/dt, ...
                                                    vertcat(model.options.trend_bands.bounds), 0);
         for iband=1:length(band_indexes)
             model.constants.trend_bands(iband).name = model.options.trend_bands(iband).name;
             model.constants.trend_bands(iband).bounds = model.options.trend_bands(iband).bounds;
             model.constants.trend_bands(iband).coeff_indexes = band_indexes{iband};
         end
    else
        error(['Uknown trend type:' model.options.trend_type]);
    end
    model.constants.T = trend_mat;
    model = init_variable(model, 'trend_coeffs', {'coeff', 'channel'}, ...
                          trend_mat' * mirror_sig_bounds(model.constants.y, ...
                                                         model.constants.npb_mirror_trend_fix));
end
model.constants.Tt = model.constants.T';
assert(size(model.constants.T, 1) == model.constants.n_samples + 2 * model.constants.npb_mirror_trend_fix);

if 1
% trend signal
trend_init = unmirror_sig_bounds(model.constants.T * model.variables.trend_coeffs, model.constants.npb_mirror_trend_fix);
if isfield(simulation, 'trend')
    model = init_variable(model, 'trend', {'time', 'channel'}, ...
                          trend_init, simulation.trend);
else
    model = init_variable(model, 'trend', {'time', 'channel'}, ...
                          trend_init);
end
end

% trend_var
if isfield(simulation, 'trend_var')
    model = init_variable(model, 'trend_var', {'scalar'}, ...
                          .1,  ...
                          simulation.trend_var);
else
    model = init_variable(model, 'trend_var', {'scalar'}, ...
                          var(model.variables.trend_coeffs(:)));
end

% Init model temporary variables -> better to avoid reallocation during sampling
% allocation of partial residuals and intermediate quantities
model.tmp.part_res_resp = zeros(size(model.constants.y));
model.tmp.part_res_trend = zeros(size(model.constants.y));
model.tmp.res = zeros(size(model.constants.y));

if isfield(simulation, 'y_stim_induced')
    model.other_true_val.stim_induced = simulation.y_stim_induced;
end

if isfield(simulation, 'trend')
    model.other_true_val.trend = simulation.trend;
end

if isfield(simulation, 'noise')
    model.other_true_val.noise = simulation.noise;
end

if isfield(simulation, 'activation_flag')
    model.other_true_val.activation_flag = simulation.activation_flag;
end

end


function model = init_variable(model, var_name, var_dnames, var_init, ...
                               true_value)
                           
    assert(isfield(model.variables, var_name));
    if nargin >= 5
        model.var_true_val.(var_name) = true_value;
    end
    if model.variables_to_estimate.(var_name) || nargin < 5
        model.variables.(var_name) = var_init;
    else
        model.variables.(var_name) = model.var_true_val.(var_name);
    end
    model.var_dnames.(var_name) = var_dnames;
end