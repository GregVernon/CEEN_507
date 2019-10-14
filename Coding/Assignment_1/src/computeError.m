function err = computeError(exactSolution, feResults, method)
ELEM = feResults.Elements;
nELEM = length(ELEM);

if method == "Exact"
    err = sym(0);
    for e = 1:nELEM
        localExactSolution = exactSolution.U(ELEM(e).LocalVariate_to_GlobalVariate);
        localApproxSolution = ELEM(e).LBasisFuns' * ELEM(e).LDOF;
        elemErr = int((localExactSolution - localApproxSolution)^2 * diff(ELEM(e).LocalVariate_to_GlobalVariate) , ELEM(e).LDomain);
        err = err + sqrt(elemErr);
    end
    err = simplify(err);
elseif method == "Gauss-Quadrature"
    err = sym(0);
    [P,W] = Quadrature("Gauss-Legendre",3);
    for e = 1:nELEM
        localExactSolution = exactSolution.U(ELEM(e).LocalVariate_to_GlobalVariate);
        localExactSolution = symfun(localExactSolution, sym('x','real'));
        localApproxSolution = ELEM(e).LBasisFuns' * ELEM(e).LDOF;
        
        errFun = abs(localExactSolution - localApproxSolution)^2 * diff(ELEM(e).LocalVariate_to_GlobalVariate);
        elemErr = sym(0);
        for ii = 1:length(P)
            elemErr = elemErr + W(ii) * errFun(P(ii));
        end
        err = err + sqrt(elemErr);
    end
    err = simplify(err);
end
end