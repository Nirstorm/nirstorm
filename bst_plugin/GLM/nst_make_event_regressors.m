function events_mat = nst_make_event_regressors(events, filter, time)
% NST_MAKE_EVENT_REGRESSORS build event design matrix from given events and filter
% 
%   EVENTS_MAT = NST_MAKE_EVENT_REGRESSORS(EVENTS, FILTER, TIME
%
%      EVENTS (struct array): 
%          events as brainstorm structure (see DB_TEMPLATE('event'))
%      FILTER (1D array of double): column vector 
%          filter for convolution
%      TIME (array of float): reference time axis
%
%      Output:
%
%      EVENTS_MAT (matrix of double): 
%          Event-induced design matrix built by convolving the events with
%          the given filter.
%          Size is nb_samples x nb_event_groups 
%          
% See also DB_TEMPLATE('event'), NST_MAKE_EVENT_TOEPLITZ_MTX 

assert(size(filter, 2) == 1); % column vector

events_tpz = nst_make_event_toeplitz_mtx(events, time, size(filter, 1));
regressors = cellfun(@(m) m*filter, events_tpz, 'UniformOutput', false);
events_mat = horzcat(regressors{:});
end