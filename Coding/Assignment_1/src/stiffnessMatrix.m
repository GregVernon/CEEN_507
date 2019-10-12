% Build local matrices for each element

for e = 1:nELEM
    for A = 1:(bFun.degree + 1)
        dNA = ELEM(e).LDerivBasisFun(A);
        for B = 1:bFun.degree + 1
            dNB = ELEM(e).LDerivBasisFun(B);
            k(A,B,e) = int(dNA*dNB,ELEM(e).Ldomain);
        end
    end
end

% Assign the local nodes of each elemental stiffness matrix to a global ID

for e = 1:nELEM
    for n1 = 1:(bFun.degree + 1)
           for n2 = 1:(bFun.degree + 1)
               gID1 = eCONN(n1,e);
               gID2 = eCONN(n2,e);
               K(gID1,gID2) = K(gID1,gID2) + k(n1,n2,e);
           end
    end
end

% Things to fix: define symbolic functions
        
        
       
            
    
    
    