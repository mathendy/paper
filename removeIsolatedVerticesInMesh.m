function [x, p] = removeIsolatedVerticesInMesh(x, p)

% i = reshape(unique(p), [], 1);
i = find( sparse(1, p, 1, 1, size(x,1)) );

i2j = zeros(size(x,1), 1);
i2j(i) = 1:numel(i);

x = x(i,:);
p = i2j(p);
