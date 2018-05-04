function success = nst_download(source_fn, dest_fn, download_msg)
global GlobalData;

if nargin < 3
    download_msg = sprintf('Downloading "%s" to "%s"...\n', source_fn, dest_fn);
end

bst_interactive = ~isempty(GlobalData) && isfield(GlobalData, 'Program') && ...
                  ~isempty(GlobalData.Program) && ...
                  (~isfield(GlobalData.Program, 'isServer') || ...
                   ~GlobalData.Program.isServer);

success = 1;
tstart = tic();
if ~isempty(strfind(source_fn, 'https:')) || ~isempty(strfind(source_fn, 'ftp:')) || ~bst_interactive
    tmp_download = tempname;
    fprintf([download_msg ' ....\n']);
    if ~isempty(strfind(source_fn, 'ftp:'))
        [ftp_site, source_ftp_rfn] = nst_split_ftp(source_fn);
        hftp = ftp(ftp_site);
        tmp_ftp = tempname;
        mget(hftp, source_ftp_rfn, tmp_ftp);
        copyfile(fullfile(tmp_ftp, source_ftp_rfn), tmp_download);
        close(hftp);
    else 
        try
            tmp_download = websave(tmp_download, source_fn);
        catch ME
            errMsg = ['Download data failed from:' 10 source_fn];
            bst_error(errMsg);
            if exist(tmp_download, 'file')
                delete(tmp_download);
            end
            success = 0;
            return
        end
    end
    dest_folder = fileparts(dest_fn);
    if ~exist(dest_folder, 'dir')
        mkdir(dest_folder);
    end
    copyfile(tmp_download, dest_fn);
elseif ~isempty(strfind(source_fn, 'http'))
    errMsg = gui_brainstorm('DownloadFile', source_fn, dest_fn, basename(dest_fn));
    % Error message
    if ~isempty(errMsg)
        errMsg = ['Download data failed from:' 10 source_fn];
        bst_error(errMsg);
        success = 0;
        return
    end
else
    errMsg = ['Unhandled web ressource:' 10 url];
    bst_error(errMsg);
    success = 0;
    return
end
fprintf('%s -- done in %1.1f sec.\n', download_msg, toc(tstart));
end
