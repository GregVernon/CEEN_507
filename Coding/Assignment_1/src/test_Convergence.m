clear
exactSolution = computeExactSolution(sym("x"));

nElements = [100];
for n = 1:length(nElements)
    disp("Iteration: " + num2str(n))
    xMin = 0;
    xMax = 1;
    nELEM = nElements(n);
    eDegree = 1;
    
    % Case A
    
    feSolution = main(nELEM,eDegree,2);
    err(n,1) = double(computeError(exactSolution,feSolution,"Exact"));
    err(n,2) = double(computeError(exactSolution,feSolution,"Gauss-Quadrature"));
    
end

figure
hold on
plot(nElements,err(:,1),'k')
plot(nElements,err(:,2),'b')
ax = gca;
ax.XScale = "Log";
ax.YScale = "Log";
