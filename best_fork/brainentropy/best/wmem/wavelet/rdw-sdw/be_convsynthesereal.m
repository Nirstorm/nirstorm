function out = be_convsynthesereal(in_V, in_W, H, G, n, Jcase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Project Name: EOG Correction Artifact
% Filename:     ConvSynthese.m
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Modified by:  Xavier Leturc
% Created on:   July 21st, 2008
% Revised on:   July 1st, 2014
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
   
    %sum1 = zeros(sz1(1),n+2*sz2(2)-1,sz1(3));
    %sum2 = zeros(sz1(1),n+2*sz2(2)-1,sz1(3));   
    n1 = sz1(2);
    n2 = size(in_W,2);
    p  = length(H);
    
    for j = 1:sz1(3)
        % On ajoute des donnees avant et apres la serie pour preparer les
        % donnees pour la convolution. La longueur des donnees ajoutees varie
        % selon les moments nulles de "Jcase".
      
        tempV= [in_V(:,n1+1-p:n1,j) in_V(:,:,j)];
        tempW= [in_W(:,n2,j) in_W(:,1:n2-1,j)];
        tempW=[tempW(:,:,j) tempW(:,1:p,j)];
            
        for i = 1:sz1(1)
            sum1(i,:,j) = filter(H(end:-1:1),1,tempV(i,:));
            sum2(i,:,j) = filter(G(end:-1:1),1,tempW(i,:));
        end
    end
 
    out = sum1(:, p+1:n1+p,:)...
        + sum2(:, p:n2+p-1,:);
    
end