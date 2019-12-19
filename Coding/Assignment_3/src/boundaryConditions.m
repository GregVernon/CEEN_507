function [K, F, BC] = boundaryConditions(K,F,BC,ELEM,BSpline,method)
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
for ii = 1:6
    val = BC.(DOFS{ii}).val;
    Loc = BC.(DOFS{ii}).location;
    [BC.(DOFS{ii}).globalNodeID,~,~] = find(isAlways(subs(Loc,symvar(Loc), nodes)));
    BC.(DOFS{ii}).globalDOF = nodeDOFS(ii,BC.(DOFS{ii}).globalNodeID);
    for e = 1:nELEM
        bFun = ELEM(e).LbFun;
        % Apply BC
        [isInElement,idx] = ismember(BC.(DOFS{ii}).globalNodeID,eCONN(:,e));
        if isInElement == true
            dNG = C{e}*ELEM(e).LDerivBasisFuns;
            dNG = formula(dNG);
            dNG = dNG(idx);
        else
            dNG = sym(0);
        end
        
        dNe = formula(ELEM(e).LDerivBasisFuns);
        JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
        nLocalNodes = length(ELEM(e).LNodes);
        
        for A = 1:nLocalNodes
            dNA = C{e}*dNe;
            dNA = dNA(A);
            
            if method == "Exact"
                gTerm = val * int(dNA*dNG * JAC, sym("xi"), ELEM(e).LDomain);
            elseif method == "GaussQuadrature"
                dNA = symfun(dNA,sym("xi"));
                dNG = symfun(dNG,sym("xi"));
                integrand = dNA * dNG * formula(JAC);
                gTerm = val * numericalQuadrature(integrand,GQ);
            end
            gTerm = formula(gTerm);
            globalDOF = nodeDOFS(ii,eCONN(A,e));%BC.(DOFS{ii}).globalNodeID);
            F(globalDOF) = F(globalDOF) - sum(gTerm); % VERIFY
        end
    end
end
    


%% Neumann BCs ( u'(x), v'(x), w'(x) )
for ii = 4:6
    val = BC.(DOFS{ii}).val;
    Loc = BC.(DOFS{ii}).location;
    BC.(DOFS{ii}).globalNodeID = find(isAlways(subs(Loc,symvar(Loc), nodes)));
    BC.(DOFS{ii}).globalDOF = nodeDOFS(ii,BC.(DOFS{ii}).globalNodeID);
    
    for e = 1:nELEM
        Ne = C{e}*formula(ELEM(e).LBasisFuns);
        % Apply h
        [isInElement,idx] = ismember(BC.(DOFS{ii}).globalNodeID,eCONN(:,e));
        if isInElement == true
            nLocalNodes = length(ELEM(e).LNodes);
            for n = 1:nLocalNodes
                if n == idx
                    Nh = Ne(n);
                    Nh = symfun(Nh,symvar(Nh)); % Turn into symbolic function
                    hTerm = Nh(ELEM(e).LNodes(n)) * val;
                    globalDOF = nodeDOFS(ii,eCONN(A,e));%BC.(DOFS{ii}).globalNodeID);
                    F(globalDOF) = F(globalDOF) - hTerm;
                end
            end
        end
    end
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