function d = nst_pdist(x,y)
% D = nst_pdist( X,Y) returns the distance between each pair of 
% observations in X and Y using the euclidean distance
% Equivalent to pdist2(x,y); Might be slower
% see also pdist2

    d = zeros(size(x,1), size(y,1));
    for i = 1:size(x,1)
        for j = 1:size(y,1)
            d(i,j) = sqrt(sum( (x(i,:) - y(j,:)).^2));
        end
    end
end
