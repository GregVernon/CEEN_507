function [K, F, BC] = boundaryConditions(K,F,BC,ELEM,BSpline)
% Extract information from B-Spline
eCONN = BSpline.elementConnectivity;
C = BSpline.decomposition.localExtractionOperator;
nodes = BSpline.nodes;
nELEM = length(C);

% Create matrix of d variables
nNodes = length(nodes);
d = sym('d', [nNodes 1]);

% Define Dirichlet (traction) and Neumann (natural) BCs
%% u(x) Boundary condition (g)
% Extract info for "g" (solution BC)
g = BC.U.val;
gLoc = BC.U.location;
% gNodeID -> Global Node ID
BC.U.gNodeID = find(isAlways(subs(gLoc,lhs(gLoc), nodes)));


%% u'(x) Boundary Condition (h)
% Extract info for "h" (derivative of solution BC)
h = BC.dU.val;
hLoc = BC.dU.location;
% gNodeID -> Global Node ID
BC.dU.gNodeID = find(isAlways(subs(hLoc,lhs(hLoc), nodes)));


%% u''(x) Boundary Condition (m)
% Extract info for "m" (2nd derivative of solution BC)
m = BC.d2U.val;
mLoc = BC.d2U.location;
% gNodeID -> Global Node ID
BC.d2U.gNodeID = find(isAlways(subs(mLoc,lhs(mLoc), nodes)));

% Subtract the m terms from the F matrix
for e = 1:nELEM
    Ne = C{e}*formula(ELEM(e).LDerivBasisFuns);
    JAC = formula(ELEM(e).Jacobian_Global_to_LocalVariate);
    % Apply m
    [isInElement,idx] = ismember(BC.d2U.gNodeID,eCONN(:,e));
    if isInElement == true
        nLocalNodes = length(ELEM(e).LNodes);
        for n = 1:nLocalNodes
            if n == idx 
                Nh = Ne(n);
                Nh = symfun(Nh,symvar(Nh)); % Turn into symbolic function
                mTerm = Nh(ELEM(e).LNodes(idx)) * m * JAC;
                F(BC.d2U.gNodeID) = F(BC.d2U.gNodeID) + mTerm;
            elseif n == idx-1
                Nh = Ne(n);
                Nh = symfun(Nh,symvar(Nh)); % Turn into symbolic function
                mTerm = Nh(ELEM(e).LNodes(idx)) * m * JAC;
                F(BC.d2U.gNodeID-1) = F(BC.d2U.gNodeID-1) + mTerm;
            end
        end
    end
end

%% u'''(x) Boundary Condition (q)
% Extract info for "q" (3rd derivative of solution BC)
q = BC.d3U.val;
qLoc = BC.d3U.location;
% gNodeID -> Global Node ID
BC.d3U.gNodeID = find(isAlways(subs(qLoc,lhs(qLoc), nodes)));

% Subtract the m terms from the F matrix
for e = 1:nELEM
    Ne = C{e}*formula(ELEM(e).LBasisFuns);
    JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
    % Apply q
    [isInElement,idx] = ismember(BC.d3U.gNodeID,eCONN(:,e));
    if isInElement == true
        nLocalNodes = length(ELEM(e).LNodes);
        for n = 1:nLocalNodes
            if n == idx
                Nh = Ne(n);
                Nh = symfun(Nh,symvar(Nh)); % Turn into symbolic function
                qTerm = Nh(ELEM(e).LNodes(n)) * q;
                F(BC.d3U.gNodeID) = F(BC.d3U.gNodeID) - qTerm;
            end
        end
    end
end
%% Update linear system structure
% Remove the gNode terms from the linear system
% Stiffness Matrix
K([BC.U.gNodeID BC.U.gNodeID+1],:) = [];
K(:,[BC.U.gNodeID BC.U.gNodeID+1]) = [];

% Force Vector
F([BC.U.gNodeID BC.U.gNodeID+1]) = [];

% Solution / Unknowns Vector
d([BC.U.gNodeID BC.U.gNodeID+1]) = [];

% Define the Kd=F relationship
K * d == F;
end