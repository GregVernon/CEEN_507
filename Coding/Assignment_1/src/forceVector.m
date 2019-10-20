function [F,f] = forceVector(ELEM,eCONN,fun)
nELEM = length(ELEM);
% Build local matrices for each element
bFun = ELEM(1).LbFun;
nLocalNodes = bFun.degree + 1;
% Exact integration method
f = sym(zeros(nLocalNodes,nELEM));
for e = 1:nELEM
    Ne = formula(ELEM(e).LBasisFuns);
    JAC = ELEM(e).Jacobian_Local_to_GlobalVariate;
    for A = 1:nLocalNodes  
        NA = Ne(A);
        NA = symfun(NA,symvar(NA));
        f(A,e) = int(NA * fun(ELEM(e).LocalVariate_to_GlobalVariate)*JAC, ELEM(e).LDomain);
    end
end

nGlobalNodes = max(max(eCONN));
F = sym(zeros(nGlobalNodes, 1));
for e = 1:nELEM
    for n = 1:nLocalNodes
        gID1 = eCONN(n,e);
        F(gID1) = F(gID1) + f(n,e); %input local values into global force vector
    end
end
