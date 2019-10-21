function result = numericalQuadrature(integrand)

assert(strcmpi(class(integrand),'symfun'))

nPoints = ceil((polynomialDegree(integrand)+1)/2);
[P,W] = Quadrature("Gauss-Legendre",nPoints);

result = sym(0);
for ii = 1:length(P)
    result = result + W(ii) * integrand(P(ii));
end
result = simplify(result);
end