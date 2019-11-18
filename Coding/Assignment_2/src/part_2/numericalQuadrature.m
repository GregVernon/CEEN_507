function result = numericalQuadrature(integrand, GQ)

assert(strcmpi(class(integrand),'symfun'))

try 
    nPoints = ceil((polynomialDegree(integrand)+1)/2);
catch
    % Non-polynomial integrand
    nPoints = 3;  % Default
end
gqRule = nPoints == [GQ(:).nPoints];
P = GQ(gqRule).P;
W = GQ(gqRule).W;

result = sym(0);
for ii = 1:length(P)
    result = result + W(ii) * integrand(P(ii));
end
result = simplify(result);
end