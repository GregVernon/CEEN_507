function [F,f] = stiffnessMatrix(ELEM,eCONN,fun)
nELEM = length(ELEM);
% Build local matrices for each element
bFun = ELEM(1).LbFun;
nLocalNodes = bFun.degree + 1;
% Exact integration method
f = sym(zeros(nLocalNodes,nELEM))
DNB = ELEM(nLocalNodes).LDerivBasisFuns(nLocalNodes)
for e = 1:nELEM
    % JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
    for A = 1:nLocalNodes
        NA = ELEM(e).LBasisFuns(A);
        DNA = ELEM(e).LDerivBasisFuns(A);
        f(A,e) = int(fun*NA,ELEM(e).LDomain)+NA(0)*h - int(DNA*DNB,ELEM(e).LDomain)*g;
        %Figure out h and g
        end
    end
end

nGlobalNodes = max(max(eCONN));
F = sym(zeros(nGlobalNodes, 1));
for e = 1:nELEM
    JAC = ELEM(e).Jacobian_Global_to_LocalVariate;
    for n1 = 1:nLocalNodes
         gID1 = eCONN(n1,e);
            F(gID1) = F(gID1) + f(n1,e)*JAC; %input local values into global force vector
    end
end
