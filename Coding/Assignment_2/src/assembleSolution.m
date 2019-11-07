function [U,ELEM] = assembleSolution(ELEM,eCONN,b,C,solutionVector)
x = sym("x","real");
d = solutionVector;
nElem = length(ELEM);

% Insert solution vector into element global degrees of freedom
eCONN = b.elementConnectivity;
for e = 1:nElem
    nLocalNodes = length(eCONN(:,e));
    for n = 1:nLocalNodes
        ELEM(e).GDOF(n) = d(eCONN(n,e));
    end
end

% Map element global degrees of freedom into local degrees of freedom
for e = 1:nElem
    nLocalNodes = length(eCONN(:,e));
    for n = 1:nLocalNodes
        ELEM(e).LDOF(n) = ELEM(e).GDOF(n) * ELEM(e).Jacobian_Local_to_GlobalVariate;
    end
end

% Create single piecewise function for ease of plotting
U = sym(NaN);
for e = 1:nElem
    nLocalNodes = length(eCONN(:,e));
    NA = C{e}*ELEM(e).LBasisFuns;
    NA = NA(ELEM(e).GlobalVariate_to_LocalVariate);
    for n = 1:nLocalNodes
        if n == 1
            U = piecewise(ELEM(e).GDomain(1)<=x<ELEM(e).GDomain(2), ELEM(e).GDOF(n) * NA(n),U);
        else
            U = piecewise(ELEM(e).GDomain(1)<=x<ELEM(e).GDomain(2), U + ELEM(e).GDOF(n) * NA(n),U);
        end
    end
end