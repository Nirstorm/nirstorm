function [fitted_model] = spree_splhblf_fit(model, not_to_track)

if nargin < 2
    not_to_track = {};
end
fitted_model = nst_algo_itmethod(model, @gs_sampling, model.max_iterations, ...
    model.var_hist_pace, @cpt_gs_observables, ...
    model.burnin_period, model.obs_hist_pace, @(m) 0, not_to_track);

%% Do some after-work computations (like rescaling):
fitted_model = model_finalize(fitted_model);
end

function [model_upd] = gs_sampling(model)

% Required matlab toolboxes:
%   - Signal Processing Toolbox ()
%   - Statistics and Machine Learning Toolbox
%   - Curve Fitting Toolbox 

if model.iteration == 5
    tic
end

% f1stats
if model.variables_to_estimate.f1
    if ~model.options.use_true_stim_induced
        model.tmp.part_res_f1 = ...
            model.constants.y - ...
            model.variables.trend; %model.constants.T * model.variables.trend_coeffs(:, ic);
    else
        model.tmp.part_res_f1 = model.other_true_val.stim_induced;
    end
    
    if 1 % version all channels at once
        %TODO: put this in tmp to gain a little time:
        pp_post_var = zeros(model.constants.f1_nb_coeffs, model.constants.f1_nb_coeffs, model.constants.n_channels);
        pp_post_mean = zeros(model.constants.n_channels, model.constants.f1_nb_coeffs);
        for ic=1:model.constants.n_channels
            inv_post_var = model.constants.PXXP ./ model.variables.noise_var(ic) ...
                + inv(model.constants.D) ./ model.variables.f1_var(ic);
            gmcy = model.constants.PtXt * model.tmp.part_res_f1(:,ic) ./ model.variables.noise_var(ic);
            pp_post_var(:,:,ic) = inv(inv_post_var);
            [R,err] = cholcov(pp_post_var(:,:,ic));
            if err ~= 0
                pp_post_var(:,:,ic) = (pp_post_var(:,:,ic) * pp_post_var(:,:,ic).')/2;
            end
            pp_post_mean(ic,:) = inv_post_var \ gmcy;
        end
        %try
            model.variables.f1 = mvnrnd(pp_post_mean, pp_post_var)';
%         catch ME
%             if (strcmp(ME.identifier,'stats:mvnrnd:BadCovariance3DSymPos'))
%                 warning('Bad sample covariance for f1');
%                 for ic=1:model.constants.n_channels
%                     pp_post_var(:,:,ic) = (pp_post_var(:,:,ic) * pp_post_var(:,:,ic).')/2;
%                 end
%                 model.variables.f1 = mvnrnd(pp_post_mean, pp_post_var)';
%             else
%                 rethrow(ME);
%             end
%        end
    else % version channel by channels
        for ic=1:model.constants.n_channels
            inv_post_var = model.constants.PXXP ./ model.variables.noise_var(ic) ...
                + inv(model.constants.D) ./ model.variables.f1_var(ic);
            gmcy = model.constants.PtXt * model.tmp.part_res_f1(:,ic) ./ model.variables.noise_var(ic);
            post_var = inv(inv_post_var);
            post_mean = inv_post_var \ gmcy;
            % use LU decomp for inversion see Paper Daunizeau & Grova IEEE TSP
            % 2005
            try
                model.variables.f1(:, ic) = mvnrnd(post_mean, post_var);
            catch ME
                if (strcmp(ME.identifier,'stats:mvnrnd:BadCovariance2DSymPos'))
                    fprintf('Warning: bad cov for f1, ic=%d, iteration=%d\n', ...
                        ic, model.iteration);
                    if 0
                        fig_from_label('phi_partial_residuals');
                        plot(model.constants.y_time_axis, model.tmp.part_res_f1(:, ic));
                        fig_from_label('resp_post_mean');
                        plot(model.constants.response_time_axis, model.constants.P1 * post_mean);
                        fig_from_label('resp_post_var');
                        imagesc(post_var);
                    end
                    post_var_fix = (post_var * post_var.')/2;
                    post_mean_fix = post_mean;
                    model.variables.f1(:, ic) = mvnrnd(post_mean_fix, post_var_fix);
                else
                    rethrow(ME);
                end
            end
        end

    end
    [peak_val, ittp] = max(model.variables.response(:, ic));
    ttp = model.constants.response_time_axis(ittp);
    
    %         edge_tol = 1; %2*model.constants.dt; %TODO: expose as parameter
    %         if ttp <= edge_tol || ttp >= model.constants.response_time_axis(end) - edge_tol
    %             ttp = model.constants.response_time_axis(end);
    %             peak_val = 0;
    %         end
    model.variables.response_time_to_peak(ic) = ttp;
    model.variables.response_max(ic) = peak_val;
    
    [undershoot_val, ittu] = min(model.variables.response(:, ic));
    model.variables.response_min(ic) = undershoot_val;
    model.variables.response_time_to_undershoot(ic) = model.constants.response_time_axis(ittu);
    model.variables.response = model.constants.P1 * model.variables.f1;
    model.variables.response_block = model.constants.X_block * model.variables.response;
end

% f1 variance
if model.variables_to_estimate.f1_var
    for ic=1:model.constants.n_channels
        fDf = model.variables.f1(:, ic)' * ...
            inv(model.constants.D) * ...
            model.variables.f1(:, ic);
        model.variables.f1_var(ic) = ...
            1 / gamrnd((model.constants.f1_nb_coeffs-1)/2, 2 / fDf);
    end
end

% noise variance
if model.variables_to_estimate.noise_var
    if 0 % -> channel-wise (slow)
        for ic=1:model.constants.n_channels
            if ~model.options.use_true_noise
                model.tmp.res = model.constants.y(:,ic) - ...
                    model.constants.XP * model.variables.f1(:, ic) - ...
                    model.variables.trend(:, ic); %model.constants.T * model.variables.trend_coeffs(:,ic);
            else
                model.tmp.res = model.other_true_val.noise(:, ic);
            end
            % assert(almost_equal(model.tmp.res, model.true_val.noise_comp(:, ic)));
            model.variables.noise_var(ic) = ...
                1 ./ gamrnd((model.constants.n_samples - 1) / 2, ...
                2 / (model.tmp.res' * model.tmp.res));
        end
    else % all channels at once -> faster
        if ~model.options.use_true_noise
            model.tmp.res = model.constants.y - ...
                model.constants.XP * model.variables.f1 - ...
                model.variables.trend; %model.constants.T * model.variables.trend_coeffs(:,ic);
        else
            model.tmp.res = model.other_true_val.noise;
        end
        model.variables.noise_var = ...
            1 ./ gamrnd((model.constants.n_samples - 1) / 2, ...
            2 ./ sum(model.tmp.res .* model.tmp.res), 1, model.constants.n_channels);
    end
end

% trend coeffs
if model.variables_to_estimate.trend_coeffs
    if 0 % -> channel wise (slow)
        if ~model.options.use_true_trend
            model.tmp.part_res_trend(:, ic) = ...
                model.constants.y(:,ic) - ...
                model.constants.XP * model.variables.f1(:, ic);
        else
            model.tmp.part_res_trend(:, ic) = model.other_true_val.trend;
        end
        post_var = 1 / (1 / model.variables.trend_var + ...
            1 / model.variables.noise_var(ic));
        post_mean = post_var .* ...
            model.constants.Tt * mirror_sig_bounds(model.tmp.part_res_trend(:, ic), model.constants.npb_mirror_trend_fix) ./ ...
            model.variables.noise_var(ic);
        model.variables.trend_coeffs(:, ic) = ...
            randn(size(model.constants.T, 2),1) * sqrt(post_var) + ...
            post_mean;
        model.variables.trend(:, ic) = ...
            unmirror_sig_bounds(model.constants.T * model.variables.trend_coeffs(:, ic), model.constants.npb_mirror_trend_fix);
    else % All channels at once (faster)
        if ~model.options.use_true_trend
            model.tmp.part_res_trend = ...
                model.constants.y - ...
                model.constants.XP * model.variables.f1;
        else
            model.tmp.part_res_trend = model.other_true_val.trend;
        end
        post_var = 1 ./ (1 / model.variables.trend_var + ...
                         1 ./ model.variables.noise_var);
        unscaled_post_mean = model.constants.Tt * mirror_sig_bounds(model.tmp.part_res_trend, model.constants.npb_mirror_trend_fix);
        post_mean = repmat(post_var ./ model.variables.noise_var, size(unscaled_post_mean,1), 1) .* unscaled_post_mean;
        model.variables.trend_coeffs = ...
            randn(size(model.constants.T, 2), model.constants.n_channels) .* sqrt(post_var) + ...
            post_mean;
        model.variables.trend = ...
            unmirror_sig_bounds(model.constants.T * model.variables.trend_coeffs, model.constants.npb_mirror_trend_fix);
    end
end

% variance of trend coeffs
if model.variables_to_estimate.trend_var
    ctc = 0.;
    for ic=1:model.constants.n_channels
        ctc = ctc + model.variables.trend_coeffs(:, ic)' * ...
            model.variables.trend_coeffs(:, ic);
    end
    model.variables.trend_var = ...
        1./gamrnd((model.constants.n_channels * size(model.constants.T, 2)-1) / 2, 2 / ctc);
end

if model.iteration == 5
    it_time = toc;
    disp(['Time for one iteration: ' num2str(it_time) ...
        ' sec. -> total should take: ' num2str(it_time * (model.max_iterations-5)) ' sec.']);
end

model_upd = model;
% plot(model.variables.response);
end

function [model] = cpt_gs_observables(model)
% Compute online posterior mean and variance for all fields
% in input *model.variables*
% TODO: push it outside this script
if ~isfield(model, 'tmp')
    model.tmp = struct();
end

field_names = fieldnames(model.variables);
for i=1:numel(field_names)
    vname = field_names{i};
    cumul_fname = ['cumul_' vname];
    cumul2_fname = ['cumul2_' vname];
    if ~isfield(model.tmp, cumul_fname)
        model.tmp.(cumul_fname) = zeros(size(model.variables.(vname)));
        model.tmp.(cumul2_fname) = zeros(size(model.variables.(vname)));
    end
    model.tmp.(cumul_fname) = model.tmp.(cumul_fname) + ...
        model.variables.(vname);
    model.tmp.(cumul2_fname) = model.tmp.(cumul2_fname) + ...
        model.variables.(vname).^2;
    
    pm_fname = [vname '_pm'];
    model.observables.(pm_fname) = model.tmp.(cumul_fname) ...
        / model.obs_iteration;
    
    pstd_fname = [vname '_pstd'];
    model.observables.(pstd_fname) = sqrt(model.tmp.(cumul2_fname) ...
        / model.obs_iteration - model.observables.(pm_fname).^2);
    
end

end

function fmodel = model_finalize(fmodel)

%% Fix samples posterior distribution which are at the edge of
% the temporal support.
extremums = {'peak', 'undershoot'};
for iext=1:length(extremums)
    sext = ['response_time_to_' extremums{iext}];
    redge_samples = fmodel.vars_history.(sext) >= ...
        fmodel.constants.response_time_axis(end) - ...
        fmodel.options.response_edge_tol;
    ledge_samples = fmodel.vars_history.(sext) <= ...
        fmodel.options.response_edge_tol;
    ic_to_fix = find(squeeze(sum(redge_samples))/length(redge_samples) >= .75 | ...
        squeeze(sum(ledge_samples))/length(ledge_samples) >= .75);
    for ichan=1:length(ic_to_fix)
        fmodel.vars_history.(sext)(ledge_samples(:, 1, ic_to_fix(ichan)), 1, ic_to_fix(ichan)) = ...
            fmodel.constants.response_time_axis(end);
        fmodel.vars_history.(sext)(redge_samples(:, 1, ic_to_fix(ichan)), 1, ic_to_fix(ichan)) = ...
            fmodel.constants.response_time_axis(end);
    end
end
%% Compute fit
fmodel.signal_fit = zeros(size(fmodel.constants.y));
fmodel.stim_induced_fit = zeros(size(fmodel.constants.y));
for ic=1:size(fmodel.signal_fit, 2)
    fmodel.stim_induced_fit(:,ic) = fmodel.constants.X * ...
        fmodel.observables.response_pm(:, ic);
    if ~isfield(fmodel.observables, 'trend_pm')
        trend = unmirror_sig_bounds(fmodel.constants.T * fmodel.observables.trend_coeffs_pm(:, ic), ...
            fmodel.constants.npb_mirror_trend_fix);
    else
        trend = fmodel.observables.trend_pm;
    end
    fmodel.signal_fit(:, ic) = fmodel.stim_induced_fit(:,ic) + trend(:,ic);
end

%% compute summary stats
bup = fmodel.burnin_period;

% a_max = max(f) [sampled random variable]
% t_a_max = dt * argmax(f) [sampled random variable]
% a_min = min(f) [sampled random variable]
% t_a_min = dt * argmin(f) [sampled random variable]
% f_pm_max = max(f^PM) [descriptor of f^PM]
% t_f_pm_max = dt * argmax(f^PM) [descriptor of f^PM]

fmodel.summary.t_a_max_pm = fmodel.observables.response_time_to_peak_pm; % \hat{\upsilon}+_pm
fmodel.summary.t_a_min_pm = fmodel.observables.response_time_to_undershoot_pm; % \hat{\upsilon}-_pm
fmodel.summary.a_max_pm = fmodel.observables.response_max_pm; % \hat{\ab}_max
fmodel.summary.a_min_pm = fmodel.observables.response_min_pm; % \hat{\ab}_min
a_max_smpls = fmodel.vars_history.response_max(bup:end, :); % p(\ab_max | y)
a_min_smpls = fmodel.vars_history.response_min(bup:end, :); % p(\ab_min | y)
t_a_max_smpls = fmodel.vars_history.response_time_to_peak(bup:end, 1, :);
t_a_min_smpls = fmodel.vars_history.response_time_to_undershoot(bup:end, :);

f_smpls = fmodel.vars_history.response(bup:end, :, :);

pk_th = fmodel.summary.pk_threshold;
alpha = 0.05; %TODO: expose alpha
percentiles = [alpha/2 1-alpha/2] .* 100;

for ichan=1:fmodel.constants.n_channels
    % Compute CI and post probas
    fmodel.summary.prob_a_max(ichan) = sum(a_max_smpls(:, ichan)  > pk_th(ichan)) / size(a_max_smpls, 1);
    fmodel.summary.t_a_max_IC(ichan,:) = prctile(t_a_max_smpls(:, ichan), percentiles);
    
    fmodel.summary.prob_a_min(ichan) = sum(a_min_smpls(:, ichan)  < -pk_th(ichan)) / size(a_min_smpls, 1);
    fmodel.summary.t_a_min_IC(ichan,:) = prctile(t_a_min_smpls(:, ichan), percentiles);
    
    % Descriptors of f_pm:
    dt_f = diff(fmodel.constants.response_time_axis(1:2));
    
    [fmodel.summary.max_f_pm(ichan), ittp] = max(fmodel.observables.response_pm(:, ichan));
    fmodel.summary.it_max_f_pm(ichan) = ittp;
    fmodel.summary.t_max_f_pm(ichan) = ittp * dt_f;
    fmodel.summary.prob_max_f_pm(ichan) = sum(f_smpls(:, ittp, ichan)  > pk_th(ichan)) / size(f_smpls, 1);
    
    [fmodel.summary.min_f_pm(ichan), ittu] = min(fmodel.observables.response_pm(:, ichan));
    fmodel.summary.it_min_f_pm(ichan) = ittu;
    fmodel.summary.t_min_f_pm(ichan) = ittu * dt_f;
    fmodel.summary.prob_min_f_pm(ichan) = sum(f_smpls(:, ittu, ichan)  < -pk_th(ichan)) / size(f_smpls, 1);
end

end
