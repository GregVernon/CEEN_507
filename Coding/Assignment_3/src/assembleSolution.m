function [U,ELEM] = assembleSolution(ELEM,BSpline,BC,solutionVector)
% Extract information from BSpline
eCONN = BSpline.elementConnectivity;
C = BSpline.decomposition.localExtractionOperator;

x = sym("x","real");
d = solutionVector;
nElem = length(ELEM);

% Node->DOF mapping
DOFS = fieldnames(BC);
nDOF = length(DOFS);
nGlobalNodes = max(max(eCONN));
nodeDOFS = reshape(1:nGlobalNodes*nDOF,nDOF,nGlobalNodes);

% Insert BCs back into solution vector
d_temp = sym(zeros(max(nodeDOFS(:)),1));

BC_DOFS = cell2mat(struct2cell(BC));
BC_DOFS = [BC_DOFS.globalDOF];

d_BC = cell2mat(struct2cell(BC));
d_BC = [d_BC.val];

d_keep = true(size(d_temp));
d_keep(BC_DOFS) = false;
d_temp(d_keep) = d;
d_temp(BC_DOFS) = d_BC;

d = d_temp;
clearvars d_BC d_keep d_temp

% Insert solution vector into element global degrees of freedom
for e = 1:nElem
    nLocalNodes = length(eCONN(:,e));
    ELEM(e).GDOF = [];
    for n = 1:nLocalNodes
        ELEM(e).GDOF(:,n) = d(nodeDOFS(:,eCONN(n,e)));
    end
end

% Map element global degrees of freedom into local degrees of freedom
for e = 1:nElem
    nLocalNodes = length(eCONN(:,e));
    ELEM(e).LDOF = [];
    for n = 1:nLocalNodes
        ELEM(e).LDOF(:,n) = ELEM(e).GDOF(:,n) * ELEM(e).Jacobian_Local_to_GlobalVariate;
    end
end

% Create single piecewise function for ease of plotting
U = sym(NaN(nDOF,1));
for e = 1:nElem
    nLocalNodes = length(eCONN(:,e));
    NA = C{e}*ELEM(e).LBasisFuns;
    NA = NA(ELEM(e).GlobalVariate_to_LocalVariate);
    for n = 1:nLocalNodes
        for ii = 1:nDOF
            if n == 1 && e == nElem
                U(ii) = piecewise(ELEM(e).GDomain(1)<=x<=ELEM(e).GDomain(2), ELEM(e).GDOF(ii,n) * NA(n),U(ii));
            elseif n == 1
                U(ii) = piecewise(ELEM(e).GDomain(1)<=x<ELEM(e).GDomain(2), ELEM(e).GDOF(ii,n) * NA(n),U(ii));
            elseif e == nElem
                U(ii) = piecewise(ELEM(e).GDomain(1)<=x<=ELEM(e).GDomain(2), U(ii) + ELEM(e).GDOF(ii,n) * NA(n),U(ii));
            else
                U(ii) = piecewise(ELEM(e).GDomain(1)<=x<ELEM(e).GDomain(2), U(ii) + ELEM(e).GDOF(ii,n) * NA(n),U(ii));
            end
        end
    end
end
U = symfun(U,symvar(U));