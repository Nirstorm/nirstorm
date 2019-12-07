function [cmat, band_indexes] = build_basis_dct(nsamples, sampling_rate, freq_ranges, ortho)
% build_basis_cosine - Build a Cosine basis within specific frequency
% ranges (or bands)
%
% Synopsis
%   [cmat] = build_basis_cosine(nsamples, sampling_rate, freq_ranges)
%
% Description
%   build an orthogonal matrix whose regressors are built from cosine
%   functions spanning specified frequency bands. In practice, building
%   such matrix consists in
%   gathering sub-columns of a Discrete Cosine Transform (DCT) matrix,
%   that corresponds to the specified frequencies.
%
% Inputs ([]s are optional)
%   (scalar int) nsamples
%       number of samples -> number of lines in the basis
%   (scalar float) sampling_rate
%       sampling rate of the signal
%   (matrix) freq_ranges
%       F x 2 matrix defining the F frequency-bands (min and max freqs)
%
% Outputs ([]s are optional)
%   (matrix)
%       nsamples x ncols matrix for the orthogonal
%       cosine basis. F is the number of frequency bands as defined by the
%       number of lines in freq_ranges.
%
% Examples: signal of 100 sec, 10 Hz sampling, physiological freq band
%           (heart beat, respiration, Meyer)
%    cmat = build_basis_cosine(1000, 10, [.2 .6; .8 1.5; .01 .1;]);
%
% TODO: assert that given freq ranges are not overlapping

% cache in tmp to avoid costly recomputations ...
% TODO: expose cached file name to clean properly after

global nirs_deconv_cache

cache_dir = '/tmp/nirs_dyna';
band_indexes = {};
if nargin < 4
    ortho = 0;
end

sfr = {};
for ifr=1:size(freq_ranges,1)
    sfr{ifr} = sprintf('%1.6f_%1.6f', freq_ranges(ifr, 1), freq_ranges(ifr, 2) );
end

idx_ranges = freq_ranges * nsamples / (sampling_rate/2);

if exist('nirs_deconv_cache', 'var') && isfield(nirs_deconv_cache, 'cmat') && ...
        nirs_deconv_cache.cmat.nsamples == nsamples && ...
        nirs_deconv_cache.cmat.sampling_rate == sampling_rate && ...
        nirs_deconv_cache.cmat.ortho == ortho
    cmat =  nirs_deconv_cache.cmat.cmat;
    band_indexes =  nirs_deconv_cache.cmat.band_indexes;
else
    cache_fn = fullfile(cache_dir, ...
        sprintf('dct_mat_%d_%1.6f_%s_%d.mat', nsamples, ...
        sampling_rate, strjoin(sfr, '_'), ortho));
    
    warning(['Cache file used for DCT matrix has to be cleaned by user: ' ...
        cache_fn]);
    
    if ~exist(cache_fn, 'file')
        cmat = dctmtx(nsamples)';
        sub_idx = [];
        cur_idx = 1;
        for fr=1:size(freq_ranges,1) % loop over freq ranges
            
            idx_ranges(2) = min(idx_ranges(2), nsamples);
            cur_idx_range = max(1, floor(idx_ranges(fr, 1))):ceil(idx_ranges(fr, 2));
            sub_idx = [sub_idx cur_idx_range];
            band_indexes{fr} = cur_idx:(cur_idx+length(cur_idx_range)-1);
            cur_idx = cur_idx + length(cur_idx_range);
        end
        cmat = cmat(:, sub_idx);
        
        if ortho
            cmat = orth(cmat);
        end
        if ~exist(cache_dir, 'dir')
            mkdir(cache_dir);
        end
        save(cache_fn, 'cmat', 'band_indexes');
    else
        cmat_data = load(cache_fn);
        cmat = cmat_data.cmat;
        band_indexes = cmat_data.band_indexes;
    end
    nirs_deconv_cache.cmat.cmat = cmat;
    nirs_deconv_cache.cmat.band_indexes =  band_indexes;
    nirs_deconv_cache.cmat.nsamples = nsamples;
    nirs_deconv_cache.cmat.ortho = ortho;
    nirs_deconv_cache.cmat.sampling_rate = sampling_rate;
end
end
