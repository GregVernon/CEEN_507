function [P,W]=Quadrature(Family,n)
if Family=="Gauss-Legendre"
    [P,W] = computeGaussLegendre(n);
end
end

function [P,W] = computeGaussLegendre(n)
    if n == 0
        P = sym(0);
        W = sym(2);
        return
    end
    bFun = basisFunction("Legendre",n,sym('x','real'));
    
    % Find the roots of the Legendre basis function
    assert(n<10) % "solve" won't return solution for polynomials >= 10
    bRoots = solve(bFun.basis==0,'Real',true,'MaxDegree',10);
    assert(isreal(bRoots)) % Make sure all roots are real
    
    % Sort the roots from min->max.  Note all roots are real
    [~,sIDX] = sort(double(bRoots));
    P = bRoots(sIDX);
    
    % Compute weights
    bFun.basis_deriv = diff(bFun.basis);
    W = 2 ./ ((1-P.^2) .* bFun.basis_deriv(P).^2);
end