function [K, F, BC] = boundaryConditions(K, F, BC, ELEM, eCONN, nodes)
%% Incorporate Boundary Conditions

% Create matrix of d variables
nNodes = length(nodes);
d = sym('d', [nNodes 1]);

% Define Dirichlet (traction) and Neumann (natural) BCs

% Extract info for "g" (solution BC)
g = BC.U.val;
gLoc = BC.U.location;
% gNodeID -> Global Node ID
BC.U.gNodeID = find(isAlways(subs(gLoc,lhs(gLoc), nodes)));

% Subtract the g terms from the F matrix
nELEM = length(ELEM);
for e = 1:nELEM
    % Apply g
    [isInElement,idx] = ismember(BC.U.gNodeID,eCONN(:,e));
    if isInElement == true
        dNG = ELEM(e).LDerivBasisFuns(idx);
    else
        dNG = sym(0);
    end
    
    dNe = formula(ELEM(e).LBasisFuns);
    JAC = ELEM(e).Jacobian_Local_to_GlobalVariate;
    nLocalNodes = length(ELEM(e).LNodes);
    for A = 1:nLocalNodes
        dNA = dNe(A);
        gTerm = int(dNA*dNG,ELEM(e).LDomain)*g;
        F(A) = F(A) - gTerm*JAC; % VERIFY
    end
end


% Extract info for "h" (derivative of solution BC)
h = BC.dU.val;
hLoc = BC.dU.location;
% gNodeID -> Global Node ID
BC.dU.gNodeID = find(isAlways(subs(hLoc,lhs(hLoc), nodes)));


% Subtract the h terms from the F matrix
for e = 1:nELEM
    %     Jac_L2G = ELEM(e).Jacobian_Local_to_GlobalVariate;
    Ne = formula(ELEM(e).LBasisFuns);
    % Apply h
    [isInElement,idx] = ismember(BC.dU.gNodeID,eCONN(:,e));
    if isInElement == true
        NH = Ne(idx);
    else
        NH = sym(0);
    end
    JAC = ELEM(e).Jacobian_Local_to_GlobalVariate;
    hTerm = NH*h; % VERIFY
    F(BC.dU.gNodeID) = F(BC.dU.gNodeID) - hTerm*JAC;
end
% Reconstruct new d matrix with g-value inserted at correct node
% d = [d(1:gBCNode-1); g; d(gBCNode+1:1)];

% Reconstruct new F and d matrices that subtract the gNode terms
% F = F - K(:,gBCNode).*d(gBCNode);
F(BC.U.gNodeID) = [];
d(BC.U.gNodeID) = [];

% Reconstruct K matrix by removing all terms associated with g BC
K(BC.U.gNodeID,:) = [];
K(:,BC.U.gNodeID) = [];

% Define the Kd=F relationship
K * d == F;
end

