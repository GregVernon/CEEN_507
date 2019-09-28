%% Main FEM function
% Simply test to see if it runs without errors
nElems = 10;
elemDegree = 1;
loadCase = 1;
assert(main(nElems,elemDegree,loadCase))

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
assert(isequal(bFun.nodes, [-1 1]));
assert(isequal(bFun.basis(1), (1-x)/2));
assert(isequal(bFun.basis(2), (1+x)/2));

x = sym('x','real');
bFun = basisFunction("Lagrange",2,x);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,2));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(isequal(bFun.nodes, [-1 0 1]));
assert(isAlways(bFun.basis(1) == (x^2 - x)/2));
assert(isAlways(bFun.basis(2) == 1-x^2));
assert(isAlways(bFun.basis(3) == (x^2 + x)/2));

x = sym('x','real');
bFun = basisFunction("Lagrange",3,x);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,3));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(isequal(bFun.nodes, [-1 -1/3 1/3 1]));
assert(isAlways(bFun.basis(1) == -( 9/16)*x^3 + (9/16)*x^2 + ( 1/16)*x - 1/16));
assert(isAlways(bFun.basis(2) ==  (27/16)*x^3 - (9/16)*x^2 - (27/16)*x + 9/16));
assert(isAlways(bFun.basis(3) == -(27/16)*x^3 - (9/16)*x^2 + (27/16)*x + 9/16));
assert(isAlways(bFun.basis(4) ==  ( 9/16)*x^3 + (9/16)*x^2 - ( 1/16)*x - 1/16));

%% Construct Finite Elements
clear
%%%%% Test simple 1 element mesh
nElems = 1;
elemDegree = 1;
bFun = basisFunction("Lagrange",1,sym('x','real'));
[eCONN,x] = generateMesh(0,1,nElems,elemDegree);
ELEM = createElements(eCONN,x,bFun);

% Test for proper field names
for ii = 1:nElems
    assert(isfield(ELEM(ii),"GDomain"));
    assert(isfield(ELEM(ii),"GNodes"));
    assert(isfield(ELEM(ii),"GDOF"));
    assert(isfield(ELEM(ii),"GBasisFuns"));
    assert(isfield(ELEM(ii),"GInterpFun"));
    
    assert(isfield(ELEM(ii),"LDomain"));
    assert(isfield(ELEM(ii),"LNodes"));
    assert(isfield(ELEM(ii),"LDOF"));
    assert(isfield(ELEM(ii),"LBasisFuns"));
    assert(isfield(ELEM(ii),"LInterpFun"));
end

% Check Global Definition
assert(isequal(ELEM(1).GDomain,[0 1]))
assert(isequal(ELEM(1).GNodes, [0 1]))
assert(isequal(ELEM(1).GDOF,[1; 1]))
assert(all(polynomialDegree(ELEM(1).GBasisFuns) == elemDegree))
assert(all(polynomialDegree(ELEM(1).GInterpFun) == elemDegree))

% Check Local Definition
assert(isequal(ELEM(1).LDomain,[-1 1]))
assert(isequal(ELEM(1).LNodes, [-1 1]))
assert(isequal(ELEM(1).LDOF,[1; 1]))
assert(all(polynomialDegree(ELEM(1).LBasisFuns) == elemDegree))
assert(all(polynomialDegree(ELEM(1).LInterpFun) == elemDegree))
