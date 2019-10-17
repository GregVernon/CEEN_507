boundaryCondition.g;
boundaryCondition.h;
boundaryCondition.gNode;

%% Incorpoorate Boundary Conditions
function [K, F, BC] = boundaryCondition(K, F, BC, ELEM, eCONN, nodes);
% Create matrix of d variables
d = sym('d', [nNodes 1]);
% Define Dirichlet (traction) and Neumann (natural) BCs
g = boundaryCondition.g;
gNode = boundaryCondition.gNode;
% Reconstruct new d matrix with g-value inserted at correct node
d = [d(1:gNode-1); g; d(gNode+1:1)];
% Reconstruct new F and d matrices that subtracts the gNode terms
F = F - K(:,gNode).*d(gNode);
F(gNode) = [];
d(gNode) = [];
% Reconstruct K matrix by removing all terms associated with g BC
K(gNode,:) = [];
K(:,gNode) = [];
% Define the Kd=F relationship
K * d == F;
end

