function [K,k] = stiffnessMatrix(ELEM,eCONN,method)
nELEM = length(ELEM);
% Build local matrices for each element
bFun = ELEM(1).LbFun;
nLocalNodes = bFun.degree + 1;
if method == "Exact"
    % Exact integration method
    k = sym(zeros(nLocalNodes,nLocalNodes,nELEM));
    for e = 1:nELEM
        % JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
        for A = 1:nLocalNodes
            dNA = ELEM(e).LDerivBasisFuns(A);
            for B = 1:nLocalNodes
                dNB = ELEM(e).LDerivBasisFuns(B);
                k(A,B,e) = int(dNA*dNB,ELEM(e).LDomain);
            end
        end
    end
elseif method == "GaussQuadrature"
    k = sym(zeros(nLocalNodes,nLocalNodes,nELEM));
    for e = 1:nELEM
        for A = 1:nLocalNodes
            dNA = ELEM(e).LDerivBasisFuns(A);
            dNA = symfun(dNA,sym("xi"));
            for B = 1:nLocalNodes
                dNB = ELEM(e).LDerivBasisFuns(B);
                dNB = symfun(dNB,sym("xi"));
                integrand = dNA*dNB;
                k(A,B,e) = numericalQuadrature(integrand);
            end
        end
    end
end

% Assign the local nodes of each elemental stiffness matrix to a global ID
nGlobalNodes = max(max(eCONN));
K = sym(zeros(nGlobalNodes));
for e = 1:nELEM
    JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
    for n1 = 1:nLocalNodes
        for n2 = 1:nLocalNodes
            gID1 = eCONN(n1,e);
            gID2 = eCONN(n2,e);
            K(gID1,gID2) = K(gID1,gID2) + k(n1,n2,e)*JAC;
        end
    end
end
