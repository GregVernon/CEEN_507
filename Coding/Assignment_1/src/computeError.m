function err = computeError(exactSolution, feResults, method)
ELEM = feResults.Elements;
nELEM = length(ELEM);

if method == "Exact"
    err = sym(0);
    for e = 1:nELEM
%         globalExactSolution = piecewise(ELEM(e).GDomain(1) <= sym("x") <= ELEM(e).GDomain(2),exactSolution.U,NaN);
        localExactSolution = exactSolution.U(ELEM(e).LocalVariate_to_GlobalVariate) * ELEM(e).Jacobian_Local_to_GlobalVariate;
        localApproxSolution = ELEM(e).LBasisFuns' * ELEM(e).LDOF;
        elemErr = int((localExactSolution - localApproxSolution)^2 , ELEM(e).GDomain);
        err = err + sqrt(elemErr);
    end
    err = simplify(err);
elseif method == "Gauss-Quadrature"
    err = sym(0);
    [P,W] = Quadrature("Gauss-Legendre",9);
    for e = 1:nELEM
        localExactSolution = exactSolution.U(ELEM(e).LocalVariate_to_GlobalVariate) * ELEM(e).Jacobian_Local_to_GlobalVariate;
        localExactSolution = symfun(localExactSolution, sym('xi','real'));
        localApproxSolution = ELEM(e).LBasisFuns' * ELEM(e).LDOF;
        
        errFun = abs(localExactSolution - localApproxSolution)^2;
        elemErr = sym(0);
        for ii = 1:length(P)
            elemErr = elemErr + W(ii) * errFun(P(ii));
        end
        err = err + sqrt(elemErr);
    end
%     err = simplify(err);
end
end