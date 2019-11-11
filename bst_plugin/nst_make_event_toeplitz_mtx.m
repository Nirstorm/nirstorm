function events_tpz = nst_make_event_toeplitz_mtx(events, time_ref, nb_coeffs)
% NST_MAKE_EVENT_TOEPLITZ_MTX build binary toeplitz matrices from given events, 
% ready for convolution.
% 
%   EVENTS_TPZ = NST_MAKE_EVENT_TOEPLITZ_MTX(EVENTS, DT, NB_SAMPLES, NB_COEFFS)
%
%      EVENTS (struct array): 
%          events as brainstorm structure (see DB_TEMPLATE('event'))
%      TIME_REF (array of float): reference temporal axis.
%      NB_COEFFS (int): number of filter coefficients
%
%      EVENTS_TPZ (cell array of matrix of double): 
%          Each matrix in the returned cell array is a binary toeplitz matrix 
%          built from one event group, with size nb_samples x nb_coeffs 
%          
% See also DB_TEMPLATE('event'), NST_MAKE_EVENT_REGRESSORS 

    assert(isscalar(nb_coeffs));
    assert(round(nb_coeffs) == nb_coeffs);
    assert(isnumeric(time_ref));
    nb_samples = length(time_ref);
   
    events_tpz = cell(1, length(events));
    for icondition=1:length(events)
        binary_seq = zeros(nb_samples, 1);
        for ievent=1:size(events(icondition).times, 2)
            i_smpl = time_to_sample_idx(events(icondition).times(:, ievent), time_ref);
            binary_seq(i_smpl(1):i_smpl(2)) = 1;
        end    
        events_tpz{icondition} = toeplitz(binary_seq, [binary_seq(1) zeros(1, nb_coeffs-1)]);
        if size(events_tpz{icondition}, 1) > nb_samples
            warning('Truncate event regressor for %s, from %d to nb_samples=%d', ...
                    events(icondition).label,  size(events_tpz{icondition}, 1), nb_samples);
            events_tpz{icondition} = events_tpz{icondition}(1:nb_samples, :);
        end
    end
end

function samples = time_to_sample_idx(time, ref_time)
if nargin < 2
    assert(all(diff(diff(time))==0));
    ref_time = time;
end
samples = round((time - ref_time(1)) / diff(ref_time(1:2))) + 1;
end
