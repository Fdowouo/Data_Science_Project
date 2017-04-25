function Ylm = compute_Ylm(lmax, costheta, phi)
%% Computes the spherical harmonics Ylm for all l and m up to lmax, at the
%  given costheta and phi.
%  Works only for vector costheta and phi. Broadcasts where possible. Cannot
%  broadcast across theta and phi.
%
%  Inputs:
%    lmax : integer
%      Maximum value of l (excluded). Ylm is computed for all l = 0...lmax-1
%      and for all m = -l...+l
%    costheta : 1D-array
%      Cosine of colatitudes (polar angles) at which Ylm must be evaluated
%    phi : 1D-array
%      Longitudes (azimuthal angles) at which Ylm must be evaluated
%
%  Returns:
%    Ylm : 3D array with lmax^2 as the size of the first dimension, and
%          broadcasted [size(costheta), size(phi)] as the sizes of the
%          second and third dimensions
%      Spherical harmonics computed for all l and m up to lmax, and at the
%      given costheta and phi
%
%  Reference:
%    Digital Library of Mathematical Functions: http://dlmf.nist.gov/14.30

% Make costheta and phi into row vectors
costheta = reshape(costheta, 1, length(costheta));
phi = reshape(phi, 1, length(phi));

% Compute associated legendre polynomials for m = 0...l
Plm = zeros([lmax^2, length(costheta)]);
for l = 0:lmax-1
	Plm(l^2+l+1:l^2+2*l+1, :) = legendre(l, costheta);
end

% Compute normalization factors
m = [];
for l = 0:lmax-1
	m = [m; [-l:l]'];
end
l = l_to_lm([0:lmax-1]);
Nlm = sqrt( (2*l+1) / (4*pi) .* factorial(l-m) ./ factorial(l+m) );
Nlm = reshape(Nlm, [lmax^2, 1]);
m = reshape(m, [lmax^2, 1]);

% Compute Ylm for m = 0...l
Ylm = bsxfun(@times, Nlm, bsxfun(@times, exp(1j * bsxfun(@times, m, phi)), Plm));

% Compute Ylm for m = -l...-1
for l = 1:lmax-1
	Ylm(l^2+l:-1:l^2+1, :) = bsxfun(@times, (-1).^[1:l]', conj(Ylm(l^2+l+2:l^2+2*l+1, :)));
end
