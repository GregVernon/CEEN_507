function bFun = basisFunction(Family, p, variate,domain)
if Family == "Lagrange"
    bFun = lagrangeBasis(p,variate,domain);
elseif Family == "Legendre"
    bFun = legendreBasis(p,variate);
elseif Family == "Bernstein"
    bFun = bernsteinBasis(p,variate,domain);
end
end

%% LAGRANGE BASIS FUNCTIONS
% A Simple Expression for Multivariate Lagrange Interpolation
% Kamron Saniee (c) SIAM, 2007
% http://evoq-eval.siam.org/Portals/0/Publications/SIURO/Vol1_Issue1/A_Simple_Expression_for_Multivariate.pdf?ver=2018-03-30-130233-050
function bFun = lagrangeBasis(p,variate,domain)
node = linspace(domain(1), domain(2), p+1);

basis = cell(p+1,1);
for ii=1:p+1  % ii is the current nodal basis function we're building
    basis{ii}(variate)=sym(1);
    for jj=1:p+1 % jj is evaluating the product series for the current node
        if ii==jj
        else
            basis{ii}(variate) = basis{ii}*((variate-node(jj))/(node(ii)-node(jj)));
        end
    end
end


bFun.name = "Lagrange";
bFun.degree = p;
bFun.variate = variate;
bFun.domain = domain;
bFun.nodes = node;
bFun.basis = transpose([basis{:}]);
end
%% LEGENDRE BASIS FUNCTIONS
% https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas
function bFun = legendreBasis(p,variate)
n = p; %Maintain consistent documentation
x = variate;

P(x) = 1/((2^n)*factorial(n))*diff((x^2-1)^n,x,n);
P(x) = simplify(P);

bFun.name = "Legendre";
bFun.degree = n;
bFun.variate = variate;
bFun.basis = P;
end
%% BERNSTEIN BASIS FUNCTIONS
% https://en.wikipedia.org/wiki/Bernstein_polynomial
% CE507 Class Notes Set No. 6 "Higher Order Finite Elements"
function bFun = bernsteinBasis(p,variate,domain)
node = linspace(domain(1), domain(2), p+1);
xi = variate;

scaleFactor = 1/((domain(2)-domain(1))^sym(p));
for a=1:p+1
    binomialCoefficient = factorial(sym(p))/(factorial(sym(a)-1)*factorial(sym(p+1-a)));
    B(a) = scaleFactor * binomialCoefficient * (domain(2)-xi)^(sym(p-(a-1)))*(-domain(1)+xi)^(sym(a)-1);
end

bFun.name = "Bernstein";
bFun.degree = p;
bFun.variate = variate;
bFun.domain = domain;
bFun.nodes = node;
bFun.basis = symfun(transpose(B),variate);
end