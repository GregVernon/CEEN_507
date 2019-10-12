function ELEM = createElements(eCONN,nodes,local_bFun)

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
    ELEM(e).GBasisFuns(globalVariate) = global_bFun.basis;
    ELEM(e).GInterpFun(globalVariate) = global_bFun.basis' * ELEM(e).GNodes';
end

% Create Local Definitions
for e = nELEM:-1:1
    ELEM(e).LDomain = local_bFun.domain;
    ELEM(e).LNodes  = local_bFun.nodes;
    ELEM(e).LDOF    = ones(length(ELEM(e).LNodes),1);
    ELEM(e).LBasisFuns(localVariate) = local_bFun.basis;
	ELEM(e).LDerivBasisFuns = diff(ELEM(e).LBasisFuns(localVariate));
    ELEM(e).LInterpFun(localVariate) = local_bFun.basis' * ELEM(e).LNodes';
end

% Create mappings ?(x) <-> x(?)
for e = nELEM:-1:1
    xi_x = computeAffineMapping(ELEM(e).LDomain, ELEM(e).GDomain, globalVariate);
    x_xi = computeAffineMapping(ELEM(e).GDomain, ELEM(e).LDomain, localVariate);
    ELEM(e).GlobalVariate_to_LocalVariate = xi_x;
    ELEM(e).LocalVariate_to_GlobalVariate = x_xi;
end

end

function F = computeAffineMapping(TargetSpace, Preimage, Domain)
A = [1 Preimage(1); 1 Preimage(2)];
b = [TargetSpace(1); TargetSpace(2)];
c = A\b;

F(Domain) = c(1) + c(2)*Domain;
end
