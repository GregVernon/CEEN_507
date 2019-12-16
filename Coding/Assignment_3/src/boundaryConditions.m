function [K, F, BC] = boundaryConditions(K,F,BC,ELEM,BSpline)
% Extract information from B-Spline
eCONN = BSpline.elementConnectivity;
C = BSpline.decomposition.localExtractionOperator;
nodes = BSpline.nodes;
nELEM = length(C);


% Node->DOF mapping
DOFS = fieldnames(BC);
nDOF = length(DOFS);
nGlobalNodes = max(max(eCONN));
nodeDOFS = reshape(1:nGlobalNodes*nDOF,nDOF,nGlobalNodes);

% Create matrix of d variables
nGlobalNodes = length(nodes);
d = sym('d', [max(nodeDOFS(:)) 1]);

% Define Dirichlet (traction) and Neumann (natural) BCs
%% Dirchlet BCs ( u(x), v(x), w(x) )
for ii = 1:3
    val = BC.(DOFS{ii}).val;
    Loc = BC.(DOFS{ii}).location;
    BC.(DOFS{ii}).globalNodeID = find(isAlways(subs(Loc,lhs(Loc), nodes)));
    BC.(DOFS{ii}).globalDOF = nodeDOFS(ii,BC.(DOFS{ii}).globalNodeID);
end

%% Neumann BCs ( u'(x), v'(x), w'(x) )
for ii = 4:6
    val = BC.(DOFS{ii}).val;
    Loc = BC.(DOFS{ii}).location;
    BC.(DOFS{ii}).globalNodeID = find(isAlways(subs(Loc,lhs(Loc), nodes)));
    if BC.(DOFS{ii}).globalNodeID == 1
        BC.(DOFS{ii}).globalNodeID = [BC.(DOFS{ii}).globalNodeID BC.(DOFS{ii}).globalNodeID+1];
    elseif BC.(DOFS{ii}).globalNodeID == nGlobalNodes
        BC.(DOFS{ii}).globalNodeID = [BC.(DOFS{ii}).globalNodeID-1 BC.(DOFS{ii}).globalNodeID];
    else
        BC.(DOFS{ii}).globalNodeID = [BC.(DOFS{ii}).globalNodeID-1 BC.(DOFS{ii}).globalNodeID BC.(DOFS{ii}).globalNodeID+1];
    end
    
    BC.(DOFS{ii}).globalDOF = nodeDOFS(ii,BC.(DOFS{ii}).globalNodeID);
    BC.(DOFS{ii}).val(1:length(BC.(DOFS{ii}).globalDOF)) = val;
end

%% Update linear system structure
% Remove the Dirichlet terms from the linear system
globalDOF = cell2mat(struct2cell(BC));
globalDOF = [globalDOF.globalDOF];
% Stiffness Matrix
K([globalDOF],:) = [];
K(:,[globalDOF]) = [];

% Force Vector
F([globalDOF]) = [];

% Solution / Unknowns Vector
d([globalDOF]) = [];

% Define the Kd=F relationship
K * d == F;
end