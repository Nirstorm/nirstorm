function out = ConvSynthese(in_V, in_W, H, G, n, Jcase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Project Name: EOG Correction Artifact
% Filename:     ConvSynthese.m
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Modified by:  Benoit Décarie, Jr Eng.
% Created on:   July 21st, 2008
% Revised on:   August 1st, 2008
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Notes:        Cette fonction permet de faire la convolution entre les
%               donnees et les filtres reels ou complexes.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Inputs
%   in_V :      vecteur ou matrice de donnees a traiter;
%   in_W :      vecteur ou matrice de donnees a traiter;
%   H :         filtre complexe H (synthese);
%   G :         filtre complexe G (synthese);
%   n :         nombre de points total;
%   Jcase :     nombre de moment nul;
% Outputs
%   out :       somme des produits de convolution avec H et G.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%START.
    sz1 = size(in_V);    sz2 = size(H);
    if (size(sz1,2) == 2)
        sz1(3) = 1;
    end
    
    sum1 = zeros(sz1(1),n+2*sz2(2)-1,sz1(3));
    sum2 = zeros(sz1(1),n+2*sz2(2)-1,sz1(3));    
    for j = 1:sz1(3)
        % On ajoute des donnees avant et apres la serie pour preparer les
        % donnees pour la convolution. La longueur des donnees ajoutees varie
        % selon les moments nulles de "Jcase".
        

       tempV= [fliplr(in_V(:,1:Jcase,j))...
            zeros(sz1(1),1) in_V(:,:,j) fliplr(in_V(:,n-Jcase-1:n-1,j))];
        tempW= [-fliplr(in_W(:,1:Jcase,j))...
           zeros(sz1(1),1) in_W(:,:,j) -fliplr(in_W(:,n-Jcase-1:n-1,j))];
       
        for i=1:sz1(1)
            sum1(i,:,j) = conv(H,tempV(i,:));
            sum2(i,:,j) = conv(G,tempW(i,:));                        
        end
    end
    out = sum1(:, sz2(2):n+sz2(2)-1,:)...
        + sum2(:, sz2(2):n+sz2(2)-1,:);
end