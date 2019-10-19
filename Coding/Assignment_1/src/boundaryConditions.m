boundaryCondition.g;
boundaryCondition.h;
boundaryCondition.gBCNode;
boundaryCondition.gELEM;

%% Incorporate Boundary Conditions
function [K, F, BC] = boundaryCondition(K, F, BC, ELEM, eCONN, nodes);

% Create matrix of d variables
d = sym('d', [nNodes 1]);

% Define Dirichlet (traction) and Neumann (natural) BCs
g = boundaryCondition.g;
gBCNode = boundaryCondition.gBCNode;
gELEM = boundaryCondition.gELEM;
nLocalNodes = bFun.degree + 1;

% Subtract the g terms from the F matrix
for e = 1:nELEM
    % Jac_L2G = ELEM(e).Jacobian_Local_to_GlobalVariate;
    dNG = ELEM(gELEM).GbFun(G);
    for A = 1:nLocalNodes
        dNA = ELEM(e).GbFun(A);
        gTerm1 = int(dNA*dNG,ELEM(e).GDomain)*g;
        F(A) = F(A) - gTerm1;
        for B = 1:nLocalNodes
            dNB = ELEM(e).GbFun(B);
            hTerm2 = int(dNB*dNG,ELEM(e).GDomain)*g;
            F(B) = F(B) - hTerm2;
        end
    end
end

% Subtract the h terms from the F matrix
for e = 1:nELEM
    Jac_L2G = ELEM(e).Jacobian_Local_to_GlobalVariate;
   for A = 1:nLocalNodes
        NA = ELEM(e).LbFun(0)*Jac_L2G;
        hTerm1 = NA*h;
        F(A) = F(A) - hTerm1;
        for B = 1:nLocalNodes
            NB = ELEM(e).LbFun(0)*Jac_L2G;
            hTerm2 = NB*h;
            F(B) = F(B) - hTerm2;
        end
    end
end

% Reconstruct new d matrix with g-value inserted at correct node
% d = [d(1:gBCNode-1); g; d(gBCNode+1:1)];

% Reconstruct new F and d matrices that subtract the gNode terms
% F = F - K(:,gBCNode).*d(gBCNode);
F(gBCNode) = [];
d(gBCNode) = [];

% Reconstruct K matrix by removing all terms associated with g BC
K(gBCNode,:) = [];
K(:,gBCNode) = [];

% Define the Kd=F relationship
K * d == F;
end

