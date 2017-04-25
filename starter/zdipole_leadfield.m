function C = zdipole_leadfield(r0)
%% Computes the spherical harmonic coefficients of the scalp EEG signal
%  corresponding to a unit dipole on the z-axis.
%
%  Inputs:
%    r0 : float
%      Radius of dipole from the center of the brain (corresponds to depth)
%
%  Returns:
%    C : array of length LMAX
%      Spherical harmonic coefficients of scalp EEG

% Dipole magnitude
dipole_mag = 1; % Unit dipole

% Radii of spheres (in cm) in the 4-sphere model
r1 = 8;   % Brain
r2 = 8.2; % CSF
r3 = 8.7; % Skull
r4 = 9.2; % Scalp

% Conductivities of different layers (w.r.t the conductivity of the brain)
sigma1 = 1;    % Normalized conductivity of the brain
sigma12 = 0.2; % Cond. of brain divided by cond. of CSF
sigma13 = 80;  % Cond. of brain divided by cond. of skull
sigma14 = 1;   % Cond. of brain divided by cond. of scalp

% Stack radii and conductivity ratios
r12 = r1/r2; r23 = r2/r3; r34 = r3/r4;
rRatio = [r12, r23, r34];
sigma34 = sigma14/sigma13; sigma23 = sigma13/sigma12;
sigmaRatio = [sigma12, sigma23, sigma34];

LMAX = 50; % Maximum bandwidth up to which we compute

% Computing gamma(4,l)'s : intermediates
for l = 1:LMAX
	gamma(4, l) = l/(l+1);
end

% Computing gamma(s,l)'s and zeta(s,l)'s : intermediates for A and B
for s = 3:-1:1
	for l = 1:LMAX
		numer = l*rRatio(s)^l - (l+1) * gamma(s+1,l)/rRatio(s)^(l+1);
		denom = rRatio(s)^l + gamma(s+1,l)/rRatio(s)^(l+1);
		zeta(s+1,l) = numer/denom;
		numer2 = (l/(l+1))*sigmaRatio(s) - zeta(s+1,l)/(l+1);
		denom2 = sigmaRatio(s) + zeta(s+1,l)/(l+1);
		gamma(s,l) = numer2/denom2;
	end
end

A = zeros(4, LMAX);
B = zeros(4, LMAX);

% Spherical harmonic oefficients on the brain surface
% (given by solving Poisson's equation)
l = 1:LMAX;
B(1, :) = dipole_mag / (4 * pi * sigma1 * r0^2) .* l .* (r0 / r1).^(l + 1);
A(1, :) = B(1, :) / gamma(1, :);

% Compute spherical harmonic coefficients of potential for each sphere
l = 0:LMAX-1;
for s = 2:4
	term(s, :) = rRatio(s-1).^l + gamma(s, :) ./ rRatio(s-1).^(l+1);
	A(s, :) = (A(s-1, :) + B(s-1, :)) ./ term(s, :);
	B(s, :) = gamma(s, :) .* A(s, :);
end

C = A(4, :) + B(4, :); % Net spherical harmonic coefficients

% Add 0 to the beginning, for the 0th coefficient
C = [0, C];
