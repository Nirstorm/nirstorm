function be_look_for_wavelab

% Find brainentropy's root
Cpath   =   which(mfilename);
idSep   =   strfind( Cpath, filesep );
Cpath   =   Cpath( 1 : idSep(end-2)-1 );

% Check if wavelab is in root
Dirs    =   dir(Cpath);
Dirs    =   {Dirs.name};
isWL    =   any( strcmpi(Dirs, 'wavelab850') );

% Process 
if ~isWL
    raise_error     =   0;
    fprintf('\n\nBEst:\tWaveLab850 was not found in your library.\n\tIt must be downloaded ...\n');
    fprintf('\tFollow these steps to fix:\n')
    fprintf('\t1) Download this package : https://statistics.stanford.edu/~wavelab/Wavelab_850/WAVELAB850.ZIP\n')
    fprintf('\t2) Unzip it\n')
    fprintf('\t3) Copy content into the folder : %s\n', Cpath)
    fprintf('\t4) Relaunch the MEM\n\n')
    error('BEst needs the wavelab toolbox. See message in the command window');
end
    
%     filename        =   fullfile(Cpath, 'wavelab850.zip'); 
%         
%     if strfind(computer, 'GLNX')
%         
%         try
%             currd       =   pwd;
%             eval(['!wget -O ' filename ' https://statistics.stanford.edu/~wavelab/Wavelab_850/WAVELAB850.ZIP']);
%             eval(['!unzip -qq ' filename ' -d ' Cpath]) 
%             
%         catch
%             raise_error = 1;
%         end
%         
%     elseif strfind(computer, 'MACI')
%         
%         try
%             filename        =   fullfile(Cpath, 'wavelab850.zip');
%             eval(['!curl -L -o ' filename ' http://www-stat.stanford.edu/~wavelab/Wavelab_850/WAVELAB850.ZIP -silent']);
%             eval(['!unzip -qq ' filename ' -d ' Cpath])
%             
%         catch
%             raise_error = 1;
%         end
%         
%     elseif strfind(computer, 'PCWIN')
%         
%         try
%             make_psscript;        
%        
%             !powershell Set-ExecutionPolicy Unrestricted
%             eval(['!powershell -ExecutionPolicy RemoteSigned -File "get_wavelab.ps1" "http://www-stat.stanford.edu/~wavelab/Wavelab_850/WAVELAB850.ZIP" ' filename]);
%             delete('get_wavelab.ps1')
%         
%             % unzip
%             pp  =   pwd;
%             cd(Cpath)
%             unzip(filename)
%             cd(pp)
%             
%         catch
%             raise_error = 1;
%         end
%     end
%     
%     delete(filename);
%     if raise_error
%         fprintf('\n\nBEst error:\tcannot download wavelab850 automatically\n\tcheck your internet security and retry or download manually\n\thttp://www-stat.stanford.edu/~wavelab/Wavelab_850/WAVELAB850.ZIP\n\tafter downloading, unzip and move folder to brainentropy root\n')
%         error('BEst needs the wavelab toolbox. See message above');
%     else
%         fprintf('\tWavelab added to brainentropy root\n')
%         warning('OFF');
%     end
%             
% end

    
return

