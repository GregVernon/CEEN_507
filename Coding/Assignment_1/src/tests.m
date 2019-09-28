%% Mesh Generation
% Create mesh "LM" array per Hughes 1.14
[eCONN,x] = generateMesh(0,2,2,1);
assert(isequal(eCONN,[1 2; 2 3]));
assert(isequal(x,[0 1 2]));

[eCONN,x] = generateMesh(0,2,2,2);
assert(isequal(eCONN,[1 2 3; 3 4 5]'));
assert(isequal(x,[0 0.5 1 1.5 2]));

[eCONN,x] = generateMesh(0,3,3,1);
assert(isequal(eCONN,[1 2; 2 3; 3 4 ]'));
assert(isequal(x,[0 1 2 3]));

[eCONN,x] = generateMesh(0,3,3,2);
assert(isequal(eCONN,[1 2 3; 3 4 5; 5 6 7]'));
assert(isequal(x,[0 0.5 1 1.5 2 2.5 3]));

%% Lagrange Basis Functions
x = sym('x','real');
bFun = basisFunction("Lagrange",1,x);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,1));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(isequal(bFun.basis(1), (1-x)/2));
assert(isequal(bFun.basis(2), (1+x)/2));

x = sym('x','real');
bFun = basisFunction("Lagrange",2,x);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,2));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(all(isAlways(bFun.basis(1) == (x^2 - x)/2)));
assert(all(isAlways(bFun.basis(2) == 1-x^2)));
assert(all(isAlways(bFun.basis(3) == (x^2 + x)/2)));

x = sym('x','real');
bFun = basisFunction("Lagrange",3,x);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,3));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(all(isAlways(bFun.basis(1) == -( 9/16)*x^3 + (9/16)*x^2 + ( 1/16)*x - 1/16)));
assert(all(isAlways(bFun.basis(2) ==  (27/16)*x^3 - (9/16)*x^2 - (27/16)*x + 9/16)));
assert(all(isAlways(bFun.basis(3) == -(27/16)*x^3 - (9/16)*x^2 + (27/16)*x + 9/16)));
assert(all(isAlways(bFun.basis(4) ==  ( 9/16)*x^3 + (9/16)*x^2 - ( 1/16)*x - 1/16)));


