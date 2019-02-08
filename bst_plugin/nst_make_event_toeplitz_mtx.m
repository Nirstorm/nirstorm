function events_tpz = nst_make_event_toeplitz_mtx(events, nb_samples, nb_coeffs)
% NST_MAKE_EVENT_TOEPLITZ_MTX build binary toeplitz matrices from given events, 
% ready for convolution.
% 
%   EVENTS_TPZ = NST_MAKE_EVENT_TOEPLITZ_MTX(EVENTS, DT, NB_SAMPLES, NB_COEFFS)
%
%      EVENTS (struct array): 
%          events as brainstorm structure (see DB_TEMPLATE('event'))
%      NB_SAMPLES (int): total number of samples
%      NB_COEFFS (int): number of filter coefficients
%
%      EVENTS_TPZ (cell array of matrix of double): 
%          Each matrix in the returned cell array is a binary toeplitz matrix 
%          built from one event group, with size nb_samples x nb_coeffs 
%          
% See also DB_TEMPLATE('event'), NST_MAKE_EVENT_REGRESSORS 

    assert(isscalar(nb_coeffs));
    assert(round(nb_coeffs) == nb_coeffs);
    assert(isscalar(nb_samples));
    assert(round(nb_samples) == nb_samples);

    events_tpz = cell(1, length(events));
    for icondition=1:length(events)
        binary_seq = zeros(nb_samples, 1);
        for ievent=1:size(events(icondition).samples, 2)
            samples = events(icondition).samples(:, ievent);
            binary_seq(samples(1):samples(2)) = 1;
        end    
        events_tpz{icondition} = toeplitz(binary_seq, [binary_seq(1) zeros(1, nb_coeffs-1)]);
        if size(events_tpz{icondition}, 1) > nb_samples
            warning('Truncate event regressor for %s, from %d to nb_samples=%d', ...
                    events(icondition).label,  size(events_tpz{icondition}, 1), nb_samples);
            events_tpz{icondition} = events_tpz{icondition}(1:nb_samples, :);
        end
    end
end
