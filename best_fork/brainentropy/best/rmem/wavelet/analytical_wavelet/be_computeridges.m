function [OPTIONS, scalo, WData] = be_computeridges( iD, ftype, OPTIONS)
% CALCUL DES RIDGES
%
% DESCRIPTION
%   Cette fonction calcule les ridges a partir de la transformee en
%   ondelettes d'un signal donne. Les cartes de ridges ne sont pas directement
%   retournes par la fonction mais attribues comme parametres des l'objet
%
%   INPUTS:
%       -   ID      :   ridge number
%       -   ftype   :   signal type
%       -   OPTIONS :   options structure (see be_main.m)
%
%   OUTPUTS:
%       -   OPTIONS :   options structure
%       -   scalo   :   scalogram
%       -   WData   :   wavelet transforme
%   
% SOURCE
%   (Aucun)
%
% EXEMPLE
%   computeRidges(c_wavelet);
% 
%% ==============================================   
% Copyright (C) 2010 - LATIS Team
%
%  Authors: LATIS team, 2010, jan 25th
%
%% ==============================================
% License 
%
% BEst is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BEst is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BEst. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------   


tic; warning off all
        
% Check data size
[ftype, OPTIONS]   =   check_size( OPTIONS, ftype );

WData   = {};
scalo   = {};
newTime	= be_closest( OPTIONS.optional.TimeSegment([1 end]),OPTIONS.mandatory.DataTime );
                
switch(ftype)
    
    case 'data'
        % Compute wavelet transform        
        for ii = 1 : numel( OPTIONS.automatic.Modality(iD).channels )
            [WData{ii}, OPTIONS]    = be_wavelet_transform( OPTIONS.automatic.Modality(iD).data(ii,:), OPTIONS, 'full' );
        end            
        t = toc;
        fprintf('\tWavelet transform done in %4.1f seconds\n', t); 
        [scalo, WData, OPTIONS]    = be_cwt2scalo(WData, OPTIONS, t);
        
    case 'baseline'
        % Compute wavelet transform
        [WData, OPTIONS]    = be_wavelet_transform( OPTIONS.automatic.Modality(iD).baseline, OPTIONS, 'full' );
        t                   = toc;
        fprintf('\tWavelet transform done in %4.1f seconds\n', t);
        scalo               = be_cwt2scalo( WData, OPTIONS, t );
        
    case 'oversize'
        % Compute wavelet tranform and scalogram
        for ii = 1 : numel( OPTIONS.automatic.Modality(iD).channels )
            [WD, OPTIONS]           = be_wavelet_transform( OPTIONS.automatic.Modality(iD).data(ii,:), OPTIONS, 'sparse' );
            [s, WD, OPTIONS, t2]    = be_cwt2scalo( {WD}, OPTIONS, [] );
            scalo(ii)               = s;
            WData{ii}               = WD{1}(:, newTime(1):newTime(2));
        end
        
        OPTIONS.automatic.TFtime    = OPTIONS.mandatory.DataTime( newTime(1):newTime(2) );
        fprintf('\tTFplanes + scalograms computed done in %4.1f seconds\n', t2)
        
    otherwise
        error('Ridge-filter: only data filetype is supported by this process.')
         
end

% Group ridges
OPTIONS     =   be_ridge_lines( be_avgcell(scalo), iD, OPTIONS);

% Filter individual ridge maps
if strcmp(OPTIONS.ridges.method,'LillyOlhede2010')
    [scalo, OPTIONS]     	= filter_individual_ridges(WData, scalo, iD, OPTIONS);
end

        
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------- HELPER FUNCTIONS --------------------------- %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ftype, OPTIONS]  =   check_size(OPTIONS, ftype)

% default
OPTIONS.automatic.TFtime  =   OPTIONS.mandatory.DataTime;

% Get data size
nbSENS  =   numel( [OPTIONS.automatic.Modality(:).channels] );
nbTIME  =   numel( OPTIONS.mandatory.DataTime );
nbLVLs  =   OPTIONS.wavelet.nb_levels+1;

% Estimate file size
EST     =   16 * nbSENS * nbTIME * nbLVLs / 1024^2;
fprintf('\tEstimated TF data size: %.2f Mb,', EST);
if EST > 1000
    fprintf(' over ''BEst'' limit. \n\t--> we compute truncated TF planes.\n')
    OPTIONS.automatic.TFcomment     =      [OPTIONS.automatic.TFcomment ' | truncated'];
    ftype   =   'oversize';
else
    ftype   =   'data';
    fprintf(' within ''BEst'' limit.\n\t--> we compute full TF planes\n');
end

return


function [SCALO, OPTIONS]  =   filter_individual_ridges(WD, SCALO, iD, OPTIONS)

% Get spectral limits
nbF     =   numel(OPTIONS.wavelet.freqs_analyzed);
lim1    =   max(1,be_closest(OPTIONS.wavelet.freqs_analyzed, OPTIONS.ridges.frequency_range(1) )-1);
lim2    =   min(nbF,be_closest(OPTIONS.wavelet.freqs_analyzed, OPTIONS.ridges.frequency_range(2) ) + 1);
newTime	=   be_closest( OPTIONS.optional.TimeSegment([1 end]),OPTIONS.mandatory.DataTime );

% Make new WD matrix
WD2     =   cellfun(@(a) full(a(lim1:lim2,:))', WD, 'uni', 0); 
clear WD
WD2     =   cat(3, WD2{:});

% Get multivariate ridges
[ir,jr]                 =   ridgewalk( WD2, fliplr(OPTIONS.wavelet.freqs_analyzed(lim1:lim2)) ); 
clear WD2
jr(ir==lim1)            =   [];
ir(ir==lim1)            =   [];
jr(ir==lim2)            =   [];
ir(ir==lim2)            =   [];
s                       =   zeros( lim2-lim1+1, diff(newTime)+1 );
idx                     =   sub2ind( size(s), jr, ir );
s(idx(~isnan(idx)))     =   1; 

% Filter inidividual ridges
for ii =1 : numel(SCALO)
SCALO{ii}(lim1:lim2,newTime(1):newTime(2)) = SCALO{ii}(lim1:lim2,newTime(1):newTime(2)).*s;
end
SKLO2                   =   zeros(size(SCALO{1}));
SKLO2(lim1:lim2,newTime(1):newTime(2)) =   s;
[Groups]                =   be_group_ridges(SKLO2, OPTIONS.ridges.min_duration, OPTIONS.ridges.scalo_threshold);
% nbN                     =   [0; find( isnan(ir) )];
% Groups                  =   {};
% for ii = 1 : numel(nbN)-1    
%     Groups{ii}          =   idx(nbN(ii)+1:nbN(ii+1)-1);
% end
OPTIONS.automatic.Modality(iD).ridges_data   =   Groups;

return
