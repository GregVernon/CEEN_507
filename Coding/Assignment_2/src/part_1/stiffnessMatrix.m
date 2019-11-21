function [K,k] = stiffnessMatrix(ELEM,BSpline,method)
% Extract information from BSpline
eCONN = BSpline.elementConnectivity;
C = BSpline.decomposition.localExtractionOperator;

nELEM = length(ELEM);

% Build local matrices for each element
bFun = ELEM(1).LbFun;
nLocalNodes = bFun.degree + 1;
if method == "Exact"
    % Exact integration method
    k = sym(zeros(nLocalNodes,nLocalNodes,nELEM));
    for e = 1:nELEM
        BS = C{e}*ELEM(e).LBasisFuns;
        dBS = diff(BS);
        for A = 1:nLocalNodes
            dNA = formula(dBS);
            dNA = dNA(A);
            for B = 1:nLocalNodes
                dNB = formula(dBS);
                dNB = dNB(B);
                k(A,B,e) = int(dNA*dNB,ELEM(e).LDomain);
            end
        end
    end
elseif method == "GaussQuadrature"
    % Precompute Gauss Quadrature rules
    ii = 0;
    for n = 9:-1:0
        ii = ii+1;
        [P,W] = Quadrature("Gauss-Legendre",n);
        GQ(ii).P = P;
        GQ(ii).W = W;
        GQ(ii).nPoints = n;
        GQ(ii).maxDegree = 2*n - 1;
    end
    
    k = sym(zeros(nLocalNodes,nLocalNodes,nELEM));
    for e = 1:nELEM
        for A = 1:nLocalNodes
            dNA = ELEM(e).LDerivBasisFuns(A);
            dNA = symfun(dNA,sym("xi"));
            for B = 1:nLocalNodes
                dNB = ELEM(e).LDerivBasisFuns(B);
                dNB = symfun(dNB,sym("xi"));
                integrand = dNA*dNB;
                k(A,B,e) = numericalQuadrature(integrand,GQ);
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
