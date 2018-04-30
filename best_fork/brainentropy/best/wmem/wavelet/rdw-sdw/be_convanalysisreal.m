function [outV outW] = be_convanalysisreal(in_C, H, G, extend, Jcase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Project Name: EOG Correction Artifact
% Filename:     ConvAnalysis.m
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Modified by:  Xavier Leturc
% Created on:   July 21st, 2008
% Revised on:   July 1st, 2014
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Notes:        Cette fonction permet de faire la convolution entre les
%               donnees et les filtres reels ou complexes.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Input
%   in_C :      vecteur ou matrice de donnees a traiter;
%   H :         filtre complexe H (analyse);
%   G :         filtre complexe G (analyse);
%   extend :    nombre de points total;
%   Jcase :     nombre de moment nul;
% Output
%   outV :      produit de convolution avec H;
%   outW :      produit de convolution avec G.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%START.
    sz1 = size(in_C);
    sz2 = size(H);
    if (size(sz1,2) == 2)
        sz1(3) = 1;
    end
    n = sz1(2);
    p = length(H);
    
    for j = 1:sz1(3)
        % On ajoute des donnees avant et apres la serie pour preparer les
        % donnees pour la convolution. La longueur des donnees ajoutees varie
        % selon les moments nulles de "Jcase".
        temp_hi=[in_C(:,2:n,j) in_C(:,1,j)];
        temp_hi=[temp_hi(:,n+1-2*Jcase-2:n,j) temp_hi(:,:,j)];
        
        temp_lo=[in_C(:,1:end,j) in_C(:,1:2*Jcase+2,j)];

        for i = 1:sz1(1)
            outV2(i,:,j) = filter(H,1,temp_lo(i,:));
            outW2(i,:,j) = filter(G,1,temp_hi(i,:));
        end
    end
    
    outV_tmp = outV2(:, 2*Jcase+2:n+2*Jcase+1,:);
    outW_tmp = outW2(:, 2*Jcase+3:n+2*Jcase+2,:);
    outV     = outV_tmp(:,1:2:(n-1),:);
    outW     = outW_tmp(:,1:2:(n-1),:);

end