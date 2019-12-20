function B = strainDisplacementMatrix(N)
%Make B a part of the basis function structure?

% Preallocate
B = sym(zeros(6,6));
val = diff(N,symvar(N)); %compute diff once
B(2,1) = val;
B(3,2) = val;
B(1,3) = val;
B(3,4) = N;
B(4,4) = val;
B(2,5) = -N;
B(5,5) = val;
B(6,6) = val;
end
