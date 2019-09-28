function bFun = basisFunction(Family, p, variate)
if Family == "Lagrange"
    bFun = lagrangeBasis(p,variate);
end
end

%% LAGRANGE BASIS FUNCTIONS
% A Simple Expression for Multivariate Lagrange Interpolation
% Kamron Saniee (c) SIAM, 2007
% http://evoq-eval.siam.org/Portals/0/Publications/SIURO/Vol1_Issue1/A_Simple_Expression_for_Multivariate.pdf?ver=2018-03-30-130233-050
function bFun = lagrangeBasis(p,variate)
domain = [-1 1];
node = linspace(domain(1), domain(2), p+1);
basis = sym('L',[p+1,1]);

for ii=1:p+1
    basis(ii)=sym(1);
    for jj=1:p+1
        if ii==jj
        else
            basis(ii) = basis(ii)*((variate-node(jj))/(node(ii)-node(jj)));
        end
    end
end


bFun.name = "Lagrange";
bFun.degree = p;
bFun.variate = variate;
bFun.domain = domain;
bFun.nodes = node;
bFun.basis = basis;
end
