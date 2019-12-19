function feSolution = main(nElem, elemDegree, elemContinuity, f, BC, base, height, E, nu, shearCorrection)
%% Define problem domain
x = sym("x","real");
xMin = 0;
xMax = 1;

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
ELEM = createElements(BSpline, bFun, base, height, E, nu, shearCorrection);

%% Construct linear system of equations
[K,k] = stiffnessMatrix(ELEM,BSpline,"Exact");
F = forceVector(ELEM,BSpline,f,"Exact");

%% Apply boundary conditions
[K,F,BC] = boundaryConditions(K,F,BC,ELEM,BSpline);

%% Solve the system of equations
d = K\F;
% d = [d(1:(BC.U.gNodeID-1))+BC.U.val; BC.U.val; BC.U.val; d(BC.U.gNodeID:end)+BC.U.val];

%% Assemble the solution
[U,ELEM,d] = assembleSolution(ELEM,BSpline,BC,d);

%% Output results
feSolution.Elements = ELEM;
feSolution.Spline = BSpline;
feSolution.U = U;
feSolution.LinearSystem.K = K;
feSolution.LinearSystem.F = F;
feSolution.LinearSystem.d = d;
end