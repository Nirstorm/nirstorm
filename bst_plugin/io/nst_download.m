function success = nst_download(source_fn, dest_fn, download_msg, bst_interactive)
global GlobalData;

if nargin < 3
    download_msg = sprintf('Downloading "%s" to "%s"...\n', source_fn, dest_fn);
end

if nargin < 4
    bst_interactive = ~isempty(GlobalData) && isfield(GlobalData, 'Program') && ...
                      ~isempty(GlobalData.Program) && ...
                      (~isfield(GlobalData.Program, 'isServer') || ...
                       ~GlobalData.Program.isServer);
end

success = 1;
tstart = tic();

[root, basename_fn, ext] = fileparts(dest_fn);
errMsg = bst_websave(dest_fn, source_fn);

tend = toc(tstart);

% Error message
if ~isempty(errMsg)
    errMsg = ['Download data failed from:' 10 source_fn];
    bst_error(errMsg);
    success = 0;
end

fprintf('%s -- done in %1.1f sec.\n', download_msg, tend);
end
