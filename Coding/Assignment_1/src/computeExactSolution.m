function exactSolution = computeExactSolution(f,g,h)
syms u(x)
eqn = diff(u,x,2) + f == 0;
Du = diff(u,x);
cond = [u(1)==g, Du(0)==h];
U(x) = dsolve(eqn,cond);

exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = [0 1];
exactSolution.U = U;

end