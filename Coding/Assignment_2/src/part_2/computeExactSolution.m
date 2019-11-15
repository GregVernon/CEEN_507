function exactSolution = computeExactSolution(E,I,f,g,h,m,q)
syms u(x)
eqn = E*I*diff(u,x,4) - f == 0;
Du = diff(u,x);
D2u = E*I*diff(u,x,2);
D3u = E*I*diff(u,x,3);
cond = [u(0)==g, Du(0)==h, D2u(1)==m, D3u(1)==q];
U(x) = dsolve(eqn,cond);

exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = [0 1];
exactSolution.U = U;

end