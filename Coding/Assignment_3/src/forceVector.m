function [F,f] = forceVector(ELEM,BSpline,fun,method)
nELEM = length(ELEM);

% Extract information from B-Spline
eCONN = BSpline.elementConnectivity;
C = BSpline.decomposition.localExtractionOperator;

% Build local matrices for each element
bFun = ELEM(1).LbFun;
nLocalNodes = bFun.degree + 1;
if method == "Exact"
    % Exact integration method
    f = cell(nELEM,1);
    for e = 1:nELEM
        Ne = C{e}*formula(ELEM(e).LBasisFuns);
        JAC = ELEM(e).Jacobian_Local_to_GlobalVariate;
        for A = 1:nLocalNodes
            NA = Ne(A);
            NA = symfun(NA,symvar(NA));
            f{e}{A} = int(NA * fun(ELEM(e).LocalVariate_to_GlobalVariate)*JAC, ELEM(e).LDomain);
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
    f = sym(zeros(nLocalNodes,nELEM));
    for e = 1:nELEM
        Ne = C{e}*formula(ELEM(e).LBasisFuns);
        JAC = ELEM(e).Jacobian_Local_to_GlobalVariate;
        for A = 1:nLocalNodes
            NA = Ne(A);
            NA = symfun(NA,symvar(NA));
            integrand = NA*fun(ELEM(e).LocalVariate_to_GlobalVariate)*JAC;
            f(A,e) = numericalQuadrature(integrand, GQ);
        end
    end
end

nDOF = size(formula(fun),1);
nGlobalNodes = max(max(eCONN));
nodeDOFS = reshape(1:nGlobalNodes*nDOF,nDOF,nGlobalNodes);

F = sym(zeros(nGlobalNodes*nDOF, 1));
for e = 1:nELEM
    for n = 1:nLocalNodes
        gID1 = nodeDOFS(:,eCONN(n,e));
        F(gID1) = F(gID1) + f{e}{n}; %input local values into global force vector
    end
end
