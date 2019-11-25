function B = strainDisplacementMatrix(N)

% Preallocate
B = sym(zeros(6,6));

B(2,1) = diff(N,symvar(N));
% continue filling the rest of B

end