function xlm = l_to_lm(xl)
%% Convert a one-dimensional vector `xl`, which is indexed by l into a vector
%  `xlm`, which is jointly indexed by l and m.

L = length(xl);
xlm = zeros(L^2, 1);

for l = 0:L-1
	xlm(l^2+1:l^2+2*l+1) = xl(l+1);
end
