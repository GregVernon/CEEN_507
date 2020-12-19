function feSolution = main(nElem, elemDegree, elemContinuity, doInterp, f, g, h)
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
if doInterp == true
    bFun = basisFunction("Bernstein", elemDegree, sym("xi","real"), [-1 1]);
else
    bFun = basisFunction("Bernstein", elemDegree, sym("xi","real"), [0 1]);
end
ELEM = createElements(BSpline, bFun, doInterp);

%% Construct linear system of equations
K = stiffnessMatrix(ELEM,BSpline,"Exact");
F = forceVector(ELEM,BSpline,f,"Exact");

%% Apply boundary conditions
[K,F,BC] = boundaryConditions(K,F,BC,ELEM,BSpline,"Exact");

%% Solve the system of equations
d = K\F;
d = [d(1:(BC.U.gNodeID-1)); BC.U.val; d(BC.U.gNodeID:end)];

%% Assemble the solution
[U,ELEM_new] = assembleSolution(ELEM,BSpline,d);

%% Compute Eigenfunction composition
[eigvec,eigval] = compute_eigenspectrum(K);
v = sym(zeros(length(d),length(eigvec)));
V = sym(zeros(length(eigval),1));
for ii = 1:length(eigval)
    v(:,ii) = real([eigvec(1:(BC.U.gNodeID-1),ii); BC.U.val; eigvec(BC.U.gNodeID:end,ii)]);
    V(ii) = real(assembleSolution(ELEM,BSpline,v(:,ii)));
end
I = compute_eigencomponents(d, v);

%% Compute Preconditioner
[P,R,C] = equilibrate(double(K));
K_pre = R*P*double(K)*C;
F_pre = R*P*double(F);

[L_pre,U_pre] = ilu(sparse(K_pre));

K_pre = inv(L_pre) * K_pre * inv(U_pre);
F_pre = inv(L_pre) * F_pre;

[eigvec_pre,eigval_pre] = compute_eigenspectrum(K_pre);
v_pre = sym(zeros(length(d),length(eigvec)));
V_pre = sym(zeros(length(eigval_pre),1));
for ii = 1:length(eigval_pre)
    v_pre(:,ii) = real([eigvec_pre(1:(BC.U.gNodeID-1),ii); BC.U.val; eigvec_pre(BC.U.gNodeID:end,ii)]);
    V_pre(ii) = real(assembleSolution(ELEM,BSpline,v_pre(:,ii)));
end
I_pre = compute_eigencomponents(d, v_pre);


%% Output results
feSolution.Elements = ELEM_new;
feSolution.Spline = BSpline;
feSolution.U = U;
feSolution.LinearSystem.K = K;
feSolution.LinearSystem.F = F;
feSolution.LinearSystem.d = d;
feSolution.Eigen.vec = v;
feSolution.Eigen.val = eigval;
feSolution.Eigen.V = V;
feSolution.Eigen.coeff = I;
feSolution.Preconditioned.K = K_pre;
feSolution.Preconditioned.F = F_pre;
feSolution.Preconditioned.d = K_pre \ F_pre;
feSolution.Preconditioned.M1 = L_pre;
feSolution.Preconditioned.M2 = U_pre;
feSolution.Preconditioned.Eigen.vec = v_pre;
feSolution.Preconditioned.Eigen.val = eigval_pre;
feSolution.Preconditioned.Eigen.V = V_pre;
feSolution.Preconditioned.Eigen.coeff = I_pre;
end