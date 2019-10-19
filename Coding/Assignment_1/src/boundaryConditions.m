function [K,F,BC] = boundaryConditions(K,F,BC,ELEM,eCONN,nodes)

g = BC.U.val;
gLoc = BC.U.location;
gNode = find(isAlways(subs(gLoc,lhs(gLoc), nodes)));
BC.U.Node = gNode;

h = BC.dU.val;
hLoc = BC.dU.location;
hNode = find(isAlways(subs(hLoc,lhs(hLoc), nodes)));
BC.dU.Node = hNode;

nElem = length(ELEM);
for e = 1:nElem
    JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
    nLocalNodes = length(ELEM(e).LNodes);
    
    % Apply g
    [isInElement,idx] = ismember(gNode,eCONN(e));
    if isInElement == true
        dNg = ELEM(e).LDerivBasisFuns(idx);
    else
        dNg = sym(0);
    end
    
    for n = 1:nLocalNodes
        dNA = ELEM(e).LDerivBasisFuns(n);
        gID = eCONN(n,e);
        F(gID) = F(gID) - g * int(dNA * dNg,ELEM(e).LDomain) * JAC;
    end
    
    % Apply h
    [isInElement,idx] = ismember(hNode,eCONN(e));
    if isInElement == true
        gID = eCONN(idx,e);
        N = ELEM(e).LBasisFuns(rhs(hLoc));
        NA = N(idx);
        F(gID) = F(gID) + (NA * h) * JAC;
    end
end

K(gNode,:) = [];
K(:,gNode) = [];
F(gNode) = [];

end
