function D = dMatrix(b, h, E, nu, shearCorrection)
% Determine if Shear Correction factors are required
if shearCorrection == "yes"
    A1 = 5/6;
    A2 = 5/6;
elseif shearCorrection == "no"
    A1 = 1;
    A2 = 1;
end

% Calculate cross-sectional area "A"
A = b * h;

% Calculate the shear modulus "G" or "mu"
mu = E / (2 * (1 + nu));

% Calculate moments of inertia in strong and weak directions of rectangular
% cross-section "I1" and "I2"
I1 = b * h^3 / 12;
I2 = h * b^3 / 12;

% Calculate the polar moment of inertia of rectangular cross-section
J = b * h * (b^2 + h^2) / 12;

% Calculate the values of the diagonal terms of the D matrix
d = [ E*A mu*A1^3 mu*A2^3 E*I1 E*I2 mu*J ];

% Build the D matrix
D = diag(d);

