function version = nst_get_version()
%NST_GET_VERSION Return the current version of NIRSTORN
version = 'github-master';
major=0;
minor=0;
patches=0;

path = fileparts(which('nst_install'));

id = fopen(fullfile(path, 'bst_plugin','VERSION'));
str=fread(id,'*char' )';

version = str(9:end-1);

end

