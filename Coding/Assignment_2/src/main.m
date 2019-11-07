function feSolution = main(nElem, elemDegree, elemContinuity, f, g, h)
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

%% Construct the finite element Spline Space
[eCONN,nodes] = generateMesh(xMin,xMax,nElem,elemDegree);

% Create the B-Spline
splineSpace.degree = elemDegree;
splineSpace.uniqueKnotsVector = linspace(xMin,xMax,nElem+1);
splineSpace.continuityVector = elemContinuity;
b = bspline(splineSpace);
% Global Bezier Extraction Operator
M = b.bezierDecomposition.globalExtractionOperator;
% Local Bezier Extraction Operators -- access as C{e}
C = b.bezierDecomposition.localExtractionOperator;

% Create the local basis functions
bFun = basisFunction("Bernstein", elemDegree, sym("xi","real"), [-1 1]);
ELEM = createElements(eCONN,nodes,C,bFun);

%% Construct linear system of equations
[K,k] = stiffnessMatrix(ELEM,eCONN,b,C,"Exact");
F = forceVector(ELEM,eCONN,b,C,f,"Exact");

%% Apply boundary conditions
[K,F,BC] = boundaryConditions(K,F,BC,ELEM,eCONN,b,C,nodes,"Exact");

%% Solve the system of equations
d = K\F;
d = [d(1:(BC.U.gNodeID-1)); BC.U.val; d(BC.U.gNodeID:end)];

%% Assemble the solution
[U,ELEM] = assembleSolution(ELEM,eCONN,b,C,d);

%% Output results
feSolution.Elements = ELEM;
feSolution.Mesh.Connectivity = eCONN;
feSolution.Mesh.Nodes = nodes;
feSolution.U = U;
feSolution.LinearSystem.K = K;
feSolution.LinearSystem.F = F;
feSolution.LinearSystem.d = d;
end