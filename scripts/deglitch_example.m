function deglitch_example()
%% Generate signal with glitches
nb_samples = 500;
nb_glitches = 10;
nb_draws = 2;
signal = randn(nb_samples, nb_draws);
gen_spikes = [ones(1, nb_draws) ; randi([3,nb_samples-2], nb_glitches-2, nb_draws) ; zeros(1, nb_draws) + nb_samples];
for idraw=1:nb_draws
    signal(gen_spikes(:, idraw), idraw) = rand(nb_glitches, 1) * 5 + 5;
end

%% Detect glitches
std_factor_thresh = 2.5;
signal_tmp = [signal(2, :) ; signal ; signal(end-1, :)]; % mirror edges
grad = diff(signal_tmp);
abs_grad = abs(grad);
std_agrad = std(abs_grad);
glitch_canditates = [zeros(1, nb_draws) ; (abs_grad > std_factor_thresh * std_agrad) .* sign(grad)];
glitch_flags = [abs(diff(glitch_canditates))==2 ; false(1, nb_draws)] & glitch_canditates;
glitch_flags = glitch_flags(2:end-1, :); % remove mirrored edges

%% Clean glitches
signal_cleaned = signal;
for idraw=1:nb_draws
    signal_cleaned(glitch_flags(:, idraw), idraw) = repmat(mean(signal(:,idraw)), sum(glitch_flags(:, idraw)), 1);
    if glitch_flags(end, idraw)
       signal_cleaned(end, idraw) = mean(signal(nb_samples-2:nb_samples-1, idraw));
       glitch_flags(end, idraw) = 0;
    end
    if glitch_flags(1, idraw)
       signal_cleaned(1, idraw) = mean(signal(2:3, idraw));
       glitch_flags(1, idraw) = 0;
    end
    % Maybe the following can be vectorized but there was smthg strange with
    % find on 2D array
    glitch_i = find(glitch_flags(:, idraw));
    signal_cleaned(glitch_i, idraw) = (signal(glitch_i-1, idraw) + signal(glitch_i+1, idraw)) ./ 2;
end

%% Plot results
figure(); hold on;
plot(signal(:, 1), 'b', 'Linewidth', 3);
plot(find(glitch_flags(:,1)), signal(glitch_flags(:,1), 1), 'Xg', 'MarkerSize', 10);
plot(signal_cleaned(:, 1), 'r', 'Linewidth', 2);
plot(gen_spikes(:, 1), signal(gen_spikes(:, 1), 1), 'ok', 'MarkerSize', 10);

figure(); hold on;
plot(signal(:, 2), 'b', 'Linewidth', 3);
plot(find(glitch_flags(:,2)), signal(glitch_flags(:,2), 2), 'Xg', 'MarkerSize', 10);
plot(signal_cleaned(:, 2), 'r', 'Linewidth', 2);
plot(gen_spikes(:, 2), signal(gen_spikes(:, 2), 2), 'ok', 'MarkerSize', 10);
end
