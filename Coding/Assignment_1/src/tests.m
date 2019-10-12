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
bFun = basisFunction("Lagrange",1,x,[-1 1]);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,1));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(isequal(bFun.nodes, [-1 1]));
assert(isequal(bFun.basis(1), (1-x)/2));
assert(isequal(bFun.basis(2), (1+x)/2));

x = sym('x','real');
bFun = basisFunction("Lagrange",2,x,[-1 1]);
assert(bFun.name == "Lagrange");
assert(isequal(bFun.degree,2));
assert(isequal(bFun.variate,x));
assert(isequal(bFun.domain,[-1 1]));
assert(isequal(bFun.nodes, [-1 0 1]));
assert(isAlways(bFun.basis(1) == (x^2 - x)/2));
assert(isAlways(bFun.basis(2) == 1-x^2));
assert(isAlways(bFun.basis(3) == (x^2 + x)/2));

x = sym('x','real');
bFun = basisFunction("Lagrange",3,x,[-1 1]);
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
%%%%% Test linear, one-element mesh %%%%%
nElems = 1;
elemDegree = 1;
bFun = basisFunction("Lagrange",elemDegree,sym('x','real'),[-1 1]);
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

%%%%% Test linear, two-element mesh %%%%%
nElems = 2;
elemDegree = 1;
bFun = basisFunction("Lagrange",elemDegree,sym('x','real'),[-1 1]);
[eCONN,x] = generateMesh(0,2,nElems,elemDegree);
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

assert(isequal(ELEM(2).GDomain,[1 2]))
assert(isequal(ELEM(2).GNodes, [1 2]))
assert(isequal(ELEM(2).GDOF,[1; 1]))
assert(all(polynomialDegree(ELEM(1).GBasisFuns) == elemDegree))
assert(all(polynomialDegree(ELEM(1).GInterpFun) == elemDegree))

% Check Local Definition
assert(isequal(ELEM(1).LDomain,[-1 1]))
assert(isequal(ELEM(1).LNodes, [-1 1]))
assert(isequal(ELEM(1).LDOF,[1; 1]))
assert(all(polynomialDegree(ELEM(1).LBasisFuns) == elemDegree))
assert(all(polynomialDegree(ELEM(1).LInterpFun) == elemDegree))

assert(isequal(ELEM(1).LDomain,[-1 1]))
assert(isequal(ELEM(1).LNodes, [-1 1]))
assert(isequal(ELEM(1).LDOF,[1; 1]))
assert(all(polynomialDegree(ELEM(1).LBasisFuns) == elemDegree))
assert(all(polynomialDegree(ELEM(1).LInterpFun) == elemDegree))

%%%%% Test quadratic, one-element mesh %%%%%
nElems = 1;
elemDegree = 2;
bFun = basisFunction("Lagrange",elemDegree,sym('x','real'),[-1 1]);
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
assert(isequal(ELEM(1).GNodes, [0 0.5 1]))
assert(isequal(ELEM(1).GDOF, [1;1;1]))
assert(all(polynomialDegree(ELEM(1).GBasisFuns) == elemDegree))
% assert(all(polynomialDegree(ELEM(1).GInterpFun) == elemDegree))

% Check Local Definition
assert(isequal(ELEM(1).LDomain,[-1 1]))
assert(isequal(ELEM(1).LNodes, [-1 0 1]))
assert(isequal(ELEM(1).LDOF,[1; 1; 1]))
assert(all(polynomialDegree(ELEM(1).LBasisFuns) == elemDegree))
% assert(all(polynomialDegree(ELEM(1).LInterpFun) == elemDegree))

%%%%% Test quadratic, two-element mesh %%%%%
nElems = 2;
elemDegree = 2;
bFun = basisFunction("Lagrange",elemDegree,sym('x','real'),[-1 1]);
[eCONN,x] = generateMesh(0,2,nElems,elemDegree);
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
assert(isequal(ELEM(1).GNodes, [0 0.5 1]))
assert(isequal(ELEM(1).GDOF, [1;1;1]))
assert(all(polynomialDegree(ELEM(1).GBasisFuns) == elemDegree))
% assert(all(polynomialDegree(ELEM(1).GInterpFun) == elemDegree))

assert(isequal(ELEM(2).GDomain,[1 2]))
assert(isequal(ELEM(2).GNodes, [1 1.5 2]))
assert(isequal(ELEM(2).GDOF, [1;1;1]))
assert(all(polynomialDegree(ELEM(2).GBasisFuns) == elemDegree))
% assert(all(polynomialDegree(ELEM(2).GInterpFun) == elemDegree))

% Check Local Definition
assert(isequal(ELEM(1).LDomain,[-1 1]))
assert(isequal(ELEM(1).LNodes, [-1 0 1]))
assert(isequal(ELEM(1).LDOF,[1; 1; 1]))
assert(all(polynomialDegree(ELEM(1).LBasisFuns) == elemDegree))
% assert(all(polynomialDegree(ELEM(1).LInterpFun) == elemDegree))

assert(isequal(ELEM(2).LDomain,[-1 1]))
assert(isequal(ELEM(2).LNodes, [-1 0 1]))
assert(isequal(ELEM(2).LDOF,[1; 1; 1]))
assert(all(polynomialDegree(ELEM(2).LBasisFuns) == elemDegree))
% assert(all(polynomialDegree(ELEM(2).LInterpFun) == elemDegree))

%% Legendre Basis Polynomials
% Degree = 0
Family = "Legendre";
p = 0;
variate = sym("x");
bFun = basisFunction(Family,p,variate);
assert(isAlways(bFun.basis==sym(1)))

% Degree = 1
Family = "Legendre";
p = 1;
variate = sym("x");
bFun = basisFunction(Family,p,variate);
assert(isAlways(bFun.basis==variate))

% Degree = 2
Family = "Legendre";
p = 2;
variate = sym("x");
bFun = basisFunction(Family,p,variate);
assert(isAlways(bFun.basis==str2sym("(1/2) * (3*x^2 - 1)")))

% Verify orthonormality
Family = "Legendre";
variate = sym("x");
for m = 0:5
    for n = 0:5
        Pm = basisFunction(Family,m,variate);
        Pn = basisFunction(Family,n,variate);
        DefInt = int(Pm.basis*Pn.basis,[-1 1]);
        if m == n
            kronDelta = 1;
        else
            kronDelta = 0;
        end
        condition = 2 / ((2*n)+1) * kronDelta;
        assert(isAlways(abs(DefInt-condition)<=eps(10)))
    end
end

%% Gauss-Legendre Quadrature
clear
Family="Gauss-Legendre";

% One-point Quadrature
n=1;
[p,w]=Quadrature(Family,n);

P = sym(0);
W = sym(2);

assert(isAlways(p==P))
assert(isAlways(w==W))

% Two-point Quadrature
n = 2;
[p,w]=Quadrature(Family,n);

P = [str2sym("-1/sqrt(3)"); ...
     str2sym("1/sqrt(3)")];
        
W = [sym(1); ...
     sym(1)];

assert(all(isAlways(p==P)))
assert(all(isAlways(w==W)))

% Three-point Quadrature
n = 3;
[p,w]=Quadrature(Family,n);

P = [str2sym("-sqrt(3/5)"); ... 
     str2sym("0"); ...
     str2sym("sqrt(3/5)")];

W = [str2sym("5/9"); ...
     str2sym("8/9"); ...
     str2sym("5/9")];
 
assert(all(isAlways(p==P)))
assert(all(isAlways(w==W)))

% Four-point Quadrature
n = 4;
[p,w]=Quadrature(Family,n);

P = [str2sym("-sqrt(3/7+2/7*sqrt(6/5))"); ...
     str2sym("-sqrt(3/7-2/7*sqrt(6/5))"); ...
     str2sym("sqrt(3/7-2/7*sqrt(6/5))"); ...
     str2sym("sqrt(3/7+2/7*sqrt(6/5))")];
         
W = [str2sym("(18-sqrt(30))/36"); ...
     str2sym("(18+sqrt(30))/36"); ...
     str2sym("(18+sqrt(30))/36"); ...
     str2sym("(18-sqrt(30))/36")];
         
assert(all(isAlways(p==P)))
assert(all(isAlways(w==W)))

% Five-point Quadrature
n = 5;
[p,w]=Quadrature(Family,n);

P = [str2sym("-1/3*sqrt(5+2*sqrt(10/7))"); ...
     str2sym("-1/3*sqrt(5-2*sqrt(10/7))"); ...
     str2sym("0"); ...
     str2sym("1/3*sqrt(5-2*sqrt(10/7))"); ...
     str2sym("1/3*sqrt(5+2*sqrt(10/7))")];
          
W = [str2sym("(322-13*sqrt(70))/900"); ...
     str2sym("(322+13*sqrt(70))/900"); ...
     str2sym("128/225"); ...
     str2sym("(322+13*sqrt(70))/900"); ...
     str2sym("(322-13*sqrt(70))/900")];

assert(all(isAlways(p==P)))
assert(all(isAlways(w==W)))
