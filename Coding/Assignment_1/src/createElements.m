function ELEM = createElements(eCONN,nodes,bFun)

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
    ELEM(e).GBasisFuns = bFun.basis;
    ELEM(e).GInterpFun = bFun.basis .* ELEM(e).GDOF;
end

% Create Local Definitions
for e = 1:nELEM
    ELEM(e).LDomain = bFun.domain;
    ELEM(e).LNodes  = bFun.nodes;
    ELEM(e).LDOF    = ones(length(ELEM(e).LNodes),1);
    ELEM(e).LBasisFuns = bFun.basis;
    ELEM(e).LInterpFun = bFun.basis .* ELEM(e).LDOF;
end


end