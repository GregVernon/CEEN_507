function ELEM = createElements(eCONN,nodes,local_bFun)

% Get information about local basis function
Family = local_bFun.name;
degree = local_bFun.degree;
variate = sym('x','real');%local_bFun.variate;

% Preallocate the array of structures
nELEM = size(eCONN,2);
for e = nELEM: -1 : 1
    ELEM(e).GDomain = [];
    ELEM(e).GNodes = [];
    ELEM(e).GDOF = [];
    ELEM(e).GBasisFuns = [];
    ELEM(e).GInterpFun = [];
    
    ELEM(e).LDomain = [];
    ELEM(e).LNodes = [];
    ELEM(e).LDOF = [];
    ELEM(e).LBasisFuns = [];
    ELEM(e).LInterpFun = [];
end

% Create Global Definitions
for e = 1:nELEM
    ELEM(e).GDomain = [nodes(eCONN(1,e)) nodes(eCONN(end,e))];
    ELEM(e).GNodes  = nodes(eCONN(:,e));
    ELEM(e).GDOF    = ones(length(ELEM(e).GNodes),1);
    global_bFun = basisFunction(Family,degree,variate,ELEM(e).GDomain);
    ELEM(e).GBasisFuns = global_bFun.basis;
    ELEM(e).GInterpFun = global_bFun.basis' * ELEM(e).GNodes';
end

% Create Local Definitions
for e = 1:nELEM
    ELEM(e).LDomain = local_bFun.domain;
    ELEM(e).LNodes  = local_bFun.nodes;
    ELEM(e).LDOF    = ones(length(ELEM(e).LNodes),1);
    ELEM(e).LBasisFuns = local_bFun.basis;
    ELEM(e).LInterpFun = local_bFun.basis .* ELEM(e).LDOF;
end


end