function success = nst_download(source_fn, dest_fn, download_msg, bst_interactive)

    if nargin < 3
        download_msg = sprintf('Downloading "%s" to "%s"...\n', source_fn, dest_fn);
    end
    
    fprintf('%s -- ', download_msg)
    
    
    success = 1;
    tstart = tic();
    errMsg = bst_websave(dest_fn, source_fn);
    tend = toc(tstart);
    
    % Error message
    if ~isempty(errMsg)
        errMsg = ['Download data failed from:' 10 source_fn];
        bst_error(errMsg);
        success = 0;
    end
    
    fprintf('done in %1.1f sec.\n', tend);
end
