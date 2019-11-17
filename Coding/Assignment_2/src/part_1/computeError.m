function err = computeError(exactSolution, feResults, method)
ELEM = feResults.Elements;
nELEM = length(ELEM);

x = sym('x');
if method == "Exact"
    err = sym(0);
    for e = 1:nELEM
        C = feResults.Spline.decomposition.localExtractionOperator{e};
        localApproxSolution = (C*ELEM(e).LBasisFuns)' * ELEM(e).LDOF;
        localExactSolution = piecewise(ELEM(e).GDomain(1)<= x <= ELEM(e).GDomain(2), exactSolution.U);
        localExactSolution = localExactSolution(ELEM(e).LocalVariate_to_GlobalVariate) * ELEM(e).Jacobian_Local_to_GlobalVariate;
        
        elemErr = int((localExactSolution - localApproxSolution)^2 , ELEM(e).GDomain);
        err = err + elemErr;
    end
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
        err = err + elemErr;
    end
end
err = sqrt(err);
end