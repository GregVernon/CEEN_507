function feSolution = main(nElem, elemDegree, elemContinuity, f, g, h, m, q, EI)
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

BC.U.location = x == 0;
BC.U.val = g;

BC.dU.location = x == 0;
BC.dU.val = h;

BC.d2U.location = x == 1;
BC.d2U.val = m;

BC.d3U.location = x == 1;
BC.d3U.val = q;

%% Construct the finite element Spline Space
% Create the B-Spline
splineSpace.degree = elemDegree;
splineSpace.uniqueKnotsVector = linspace(xMin,xMax,nElem+1);
splineSpace.continuityVector = elemContinuity;
BSpline = bspline(splineSpace);

% Global Bezier Extraction Operator
M = BSpline.decomposition.globalExtractionOperator;

% Local Bezier Extraction Operators -- access as C{e}
C = BSpline.decomposition.localExtractionOperator;

% Create the local basis functions
bFun = basisFunction("Bernstein", elemDegree, sym("xi","real"), [-1 1]);
ELEM = createElements(BSpline,bFun,EI);

%% Construct linear system of equations
[K,k] = stiffnessMatrix(ELEM,BSpline,"GaussQuadrature");
F = forceVector(ELEM,BSpline,f,"GaussQuadrature");

%% Apply boundary conditions
[K,F,BC] = boundaryConditions(K,F,BC,ELEM,BSpline);

%% Solve the system of equations
d = K\F;
d = [d(1:(BC.U.gNodeID-1))+BC.U.val; BC.U.val; BC.U.val; d(BC.U.gNodeID:end)+BC.U.val];

%% Assemble the solution
[U,ELEM] = assembleSolution(ELEM,BSpline,d);

%% Output results
feSolution.Elements = ELEM;
feSolution.Spline = BSpline;
feSolution.U = U;
feSolution.LinearSystem.K = K;
feSolution.LinearSystem.F = F;
feSolution.LinearSystem.d = d;
end