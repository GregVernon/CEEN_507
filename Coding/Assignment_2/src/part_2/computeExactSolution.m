function exactSolution = computeExactSolution(E,I,f,g,h,m,q)
syms u(x)
eqn = diff(u,x,4) == f/(E*I);
Du = diff(u,x);
D2u = diff(u,x,2);
D3u = diff(u,x,3);
cond = [u(0)==g, Du(0)==h, D2u(1)==m/(E*I), D3u(1)==q/(E*I)];
U(x) = dsolve(eqn,cond);

exactSolution.f = f;
exactSolution.eqn = eqn;
exactSolution.cond = cond;
exactSolution.domain = [0 1];
exactSolution.U = U;

end