function feSolution = main(nElem, elemDegree, f, g, h)
%% Define problem domain
x = sym("x","real");
xMin = 0;
xMax = 1;

%% Define Load-Case
%%%%% Example inputs
% if loadCase == 1
%     q = sym(1);
%     f(x) = q;
% elseif loadCase == 2
%     f(x) = x;
% elseif loadCase == 3
%     f(x) = x^2;
% end
% 
%% Define boundary conditions
% g = sym(1); 
% h = sym(0.25);

BC.U.location = x==1;
BC.U.val = g;
BC.dU.location = x == 0;
BC.dU.val = h;

%% Construct the finite element space
[eCONN,nodes] = generateMesh(xMin,xMax,nElem,elemDegree);
bFun = basisFunction("Lagrange", elemDegree, sym("xi","real"), [-1 1]);
ELEM = createElements(eCONN,nodes,bFun);

%% Construct linear system of equations
K = stiffnessMatrix(ELEM,eCONN,"GaussQuadrature");
F = forceVector(ELEM,eCONN,f,"GaussQuadrature");

%% Apply boundary conditions
[K,F,BC] = boundaryConditions(K,F,BC,ELEM,eCONN,nodes,"GaussQuadrature");

%% Solve the system of equations
d = K\F;
d = [d(1:(BC.U.gNodeID-1)); BC.U.val; d(BC.U.gNodeID:end)];

%% Assemble the solution
[U,ELEM] = assembleSolution(ELEM,eCONN,d);

%% Output results
feSolution.Elements = ELEM;
feSolution.Mesh.Connectivity = eCONN;
feSolution.Mesh.Nodes = nodes;
feSolution.U = U;
feSolution.LinearSystem.K = K;
feSolution.LinearSystem.F = F;
feSolution.LinearSystem.d = d;
end