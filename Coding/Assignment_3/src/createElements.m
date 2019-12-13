function ELEM = createElements(BSpline,local_bFun,EI)
% Extract information from B-Spline
eCONN = BSpline.decomposition.spline.elementConnectivity;
nodes = BSpline.decomposition.spline.nodes;

% Get information about local basis function
Family = local_bFun.name;
degree = local_bFun.degree;
localVariate = local_bFun.variate;
globalVariate = sym('x','real');

% Preallocate the array of structures
nELEM = size(eCONN,2);

% Create Global Definitions
for e = nELEM:-1:1
    ELEM(e).GDomain = [nodes(eCONN(1,e)) nodes(eCONN(end,e))];
    ELEM(e).GNodes  = nodes(eCONN(:,e));
    ELEM(e).GDOF    = ones(length(ELEM(e).GNodes),1);
    global_bFun = basisFunction(Family,degree,globalVariate,ELEM(e).GDomain);
    ELEM(e).GbFun = global_bFun;
    ELEM(e).GBasisFuns(globalVariate) = global_bFun.basis;
    ELEM(e).GInterpFun(globalVariate) = global_bFun.basis' * ELEM(e).GNodes;
    ELEM(e).G_EI = EI;
end

% Create Local Definitions
for e = nELEM:-1:1
    ELEM(e).LDomain = local_bFun.domain;
    ELEM(e).LNodes  = local_bFun.nodes;
    ELEM(e).LDOF    = ones(length(ELEM(e).LNodes),1);
    ELEM(e).LbFun = local_bFun;
    ELEM(e).LBasisFuns(localVariate) = local_bFun.basis;
	ELEM(e).LDerivBasisFuns(localVariate) = diff(ELEM(e).LBasisFuns(localVariate));
    ELEM(e).LInterpFun(localVariate) = local_bFun.basis' * ELEM(e).LNodes';
    ELEM(e).L_D = strainDisplacementMatrix(symfun(ELEM(e).LBasisFuns));
end

% Create mappings ?(x) <-> x(?)
for e = nELEM:-1:1
    xi_x = computeAffineMapping(ELEM(e).LDomain, ELEM(e).GDomain, globalVariate);
    x_xi = computeAffineMapping(ELEM(e).GDomain, ELEM(e).LDomain, localVariate);
    ELEM(e).GlobalVariate_to_LocalVariate = xi_x;
    ELEM(e).LocalVariate_to_GlobalVariate = x_xi;
    ELEM(e).Jacobian_Global_to_LocalVariate = diff(xi_x);
    ELEM(e).Jacobian_Local_to_GlobalVariate = diff(x_xi);
    ELEM(e).Hessian_Global_to_LocalVariate = diff(xi_x,2);
    ELEM(e).Hessian_Local_to_GlobalVariate = diff(x_xi,2);
end

end

function F = computeAffineMapping(TargetSpace, Preimage, Domain)
A = [1 Preimage(1); 1 Preimage(2)];
b = [TargetSpace(1); TargetSpace(2)];
c = A\b;

F(Domain) = c(1) + c(2)*Domain;
end
