function [K, F, BC] = boundaryConditions(K, F, BC, ELEM, eCONN, nodes, method)
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
    bFun = ELEM(e).LbFun;
    % Apply g
    [isInElement,idx] = ismember(BC.U.gNodeID,eCONN(:,e));
    if isInElement == true
        dNG = ELEM(e).LDerivBasisFuns(idx);
    else
        dNG = sym(0);
    end
    
    dNe = formula(ELEM(e).LDerivBasisFuns);
    JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
    nLocalNodes = length(ELEM(e).LNodes);
    for A = 1:nLocalNodes
        dNA = dNe(A);
        if method == "Exact"
            gTerm = g * int(dNA*dNG * JAC, sym("xi"), ELEM(e).LDomain);
        elseif method == "GaussQuadrature"
            nPoints = ceil((bFun.degree+1)/2);
            [P,W] = Quadrature("Gauss-Legendre",nPoints);
            gTerm = sym(0);
            dNA = symfun(dNA,sym("xi"));
            dNG = symfun(dNG,sym("xi"));
            integrand = dNA * dNG * formula(JAC);
            for ii = 1:length(P)
                gTerm = gTerm + g * W(ii) * integrand(P(ii));
            end
        end
        globalNodeID = eCONN(A,e);
        F(globalNodeID) = F(globalNodeID) - gTerm; % VERIFY
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
        nLocalNodes = length(ELEM(e).LNodes);
        for n = 1:nLocalNodes
            if n == idx
                Nh = Ne(n);
                Nh = symfun(Nh,symvar(Nh)); % Turn into symbolic function
                hTerm = Nh(ELEM(e).LNodes(n)) * h;
                F(BC.dU.gNodeID) = F(BC.dU.gNodeID) - hTerm;
            end
        end
    end
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