function [version_str] = rdate_to_version(rdate)
% Return the version number based on the given release date
% Eg rdate_to_version('2013a') returns '8.1'
% 
% Raise error if given release date cannot be converted

switch(rdate)
    case 'R2006a'
        version_str = '7.2';
    case 'R2006b'
        version_str = '7.3';
    case 'R2007a'
        version_str = '7.4';
    case 'R2007b'
        version_str = '7.5';
    case 'R2008a'
        version_str = '7.6';
    case 'R2008b'
        version_str = '7.7';
    case 'R2009a'
        version_str = '7.8';
    case 'R2009b'
        version_str = '7.9';
    case 'R2009bSP1'
        version_str = '7.9.1';
    case 'R2010a'
        version_str = '7.10';
    case 'R2010b'
        version_str = '7.11';
    case 'R2010bSP1'
        version_str = '7.11.1';
    case 'R2010bSP2'
        version_str = '7.11.2';
    case 'R2011a'
        version_str = '7.12';
    case 'R2011b'
        version_str = '7.13';
    case 'R2012a'
        version_str = '7.14';
    case 'R2013a'
        version_str = '8.1';
    case 'R2013b'
        version_str = '8.2';
    case 'R2014a'
        version_str = '8.3';
    case 'R2014b'
        version_str = '8.4';
    case 'R2015a'
        version_str = '8.5';
    case 'R2015aSP1'
        version_str = '8.5';
    case 'R2015b'
        version_str = '8.6';
    case 'R2016a'
        version_str = '9.0';
    case 'R2016b'
        version_str = '9.1';
    case 'R2017a'
        version_str = '9.2';
    case 'R2017b'
        version_str = '9.3';
    case 'R2018b'
        version_str = '9.4';        
    otherwise
        throw(MException('Nirstorm:BadMatlabVersionString', ['Unknown release date: ' rdate]));
end