function exactSolution = computeExactSolution(L,E,I,k,A,G,f,conditions)
syms u(x)
eqn = diff(u,x,4) == f/(E*I) - 1/(k*A*G) * diff(f,x,2);

% conditions = [deriv_1, x_1, val_1;...
%               deriv_2, x_2, val_2;...
%               deriv_3, x_3, val_3;...
%               deriv_4, x_4, val_4]

assert(isequal(size(conditions),[4,3]))
cond = sym(zeros(1,4));
for ii = 1:size(conditions,1)
    deriv_i = conditions(ii,1);
    x_i = conditions(ii,2);
    val_i = conditions(ii,3);
    du = diff(u,x,deriv_i);
    cond(ii) = du(x_i)==val_i;
end
U(x) = dsolve(eqn,cond);

exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = [0 L];
exactSolution.U = U;

end