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
end