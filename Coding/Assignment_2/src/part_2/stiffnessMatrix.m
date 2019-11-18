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
        BS = symfun(BS,symvar(BS));
        d2BS = diff(BS,2);
        d2NS = d2BS * ELEM(e).Jacobian_Local_to_GlobalVariate^(-2);
        d2NS = formula(d2NS);
        for A = 1:nLocalNodes
            d2NA = d2NS(A);
            for B = 1:nLocalNodes
                d2NB = d2NS(B);
                k(A,B,e) = int(d2NA*d2NB*ELEM(e).G_EI,ELEM(e).LDomain);
                k(A,B,e) = k(A,B,e) * ELEM(e).Jacobian_Local_to_GlobalVariate;
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
        BS = C{e}*ELEM(e).LBasisFuns;
        BS = symfun(BS,symvar(BS));
        d2BS = diff(BS,2);
        d2NS = d2BS * ELEM(e).Jacobian_Local_to_GlobalVariate^(-2);
        d2NS = formula(d2NS);
        for A = 1:nLocalNodes
            d2NA = symfun(d2NS(A),symvar(BS));
            for B = 1:nLocalNodes
                d2NB = symfun(d2NS(B),symvar(BS));
                integrand = d2NA*d2NB*ELEM(e).G_EI;
                k(A,B,e) = numericalQuadrature(integrand,GQ) * ELEM(e).Jacobian_Local_to_GlobalVariate;
            end
        end
    end
end

% Assign the local nodes of each elemental stiffness matrix to a global ID
nGlobalNodes = max(max(eCONN));
K = sym(zeros(nGlobalNodes));
for e = 1:nELEM
    for n1 = 1:nLocalNodes
        for n2 = 1:nLocalNodes
            gID1 = eCONN(n1,e);
            gID2 = eCONN(n2,e);
            K(gID1,gID2) = K(gID1,gID2) + k(n1,n2,e);
        end
    end
end
