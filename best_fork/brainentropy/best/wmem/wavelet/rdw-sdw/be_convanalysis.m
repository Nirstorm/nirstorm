function [outV outW] = ConvAnalysis(in_C, H, G, extend, Jcase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Project Name: EOG Correction Artifact
% Filename:     ConvAnalysis.m
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Modified by:  Benoit Décarie, Jr Eng.
% Created on:   July 21st, 2008
% Revised on:   August 1st, 2008
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
    sz1 = size(in_C);    sz2 = size(H);
    if (size(sz1,2) == 2)
        sz1(3) = 1;
    end
    
    outV2 = zeros(sz1(1),extend+2*sz2(2)-1,sz1(3));
    outW2 = zeros(sz1(1),extend+2*sz2(2)-1,sz1(3));
    for j = 1:sz1(3)
        % On ajoute des donnees avant et apres la serie pour preparer les
        % donnees pour la convolution. La longueur des donnees ajoutees varie
        % selon les moments nulles de "Jcase".
        temp = [fliplr(in_C(:,1:Jcase+1,j))...
                in_C(:,1:extend,j)...
                fliplr(in_C(:,extend-Jcase:extend,j))];        
        for i = 1:sz1(1)
            outV2(i,:,j) = conv(H,temp(i,:));
            outW2(i,:,j) = conv(G,temp(i,:));
        end
    end
    outV = outV2(:, sz2(2)+1:2:extend+sz2(2)-1,:);
    outW = outW2(:, sz2(2)+1:2:extend+sz2(2)-1,:);

end