b = 0.005; %problem definition
H = 0.005; %problem definition
L = 1;
E = 1000000;


syms f(x) g h
f(x) = cantileverXDeflection(H,b,L,E); %fix for beam problem
g = sym(0);
h = sym(0);

exactSolution = computeExactSolution(f,g,h);
hold on
fplot(exactSolution.U,exactSolution.domain)

solutionType = "Exact";
Error = 0;
NumElements =  inf;
ElementDegree = NaN;

save("exact.mat",'exactSolution', 'solutionType', 'Error', 'NumElements', 'ElementDegree', 'f', 'g', 'h');
exSolution = exactSolution;
exactSolution = exSolution;


nELEM = 10.^(0:2);
eDegree = 2; %Quadratic splines
for n = 1:length(nELEM)
    iter = iter+1;
    disp("Iteration " + num2str(n) + " of " + num2str(length(nELEM)) + "; Number of Elements = " + num2str(nELEM(n)))
    
    elemContinuity = [-1 (eDegree-1)*ones(1,nELEM(n)-1) -1];
    feSolution = main(nELEM(n),eDegree,elemContinuity,f,g,h);
    %         fplot(feSolution.U,[0 1])
    err = computeError(exactSolution, feSolution, "Exact");
    
    solutionType = "Finite Element - Exact Integration";
    Error = err;
    NumElements =  nELEM(n);
    ElementDegree = eDegree;
    save("approx_"+num2str(iter)+".mat","feSolution", "solutionType", "Error", "NumElements", "ElementDegree", "f", "g", "h");
    fprintf(fOut,"%s,%.16E,%d,%d,%s,%s,%s\n",solutionType,Error,NumElements,ElementDegree,string(f),string(g),string(h));
end


nELEM = 10.^(0:2);
eDegree = 3; %Cubic splines
for n = 1:length(nELEM)
    iter = iter+1;
    disp("Iteration " + num2str(n) + " of " + num2str(length(nELEM)) + "; Number of Elements = " + num2str(nELEM(n)))
    
    elemContinuity = [-1 (eDegree-1)*ones(1,nELEM(n)-1) -1];
    feSolution = main(nELEM(n),eDegree,elemContinuity,f,g,h);
    %         fplot(feSolution.U,[0 1])
    err = computeError(exactSolution, feSolution, "Exact");
    
    solutionType = "Finite Element - Exact Integration";
    Error = err;
    NumElements =  nELEM(n);
    ElementDegree = eDegree;
    save("approx_"+num2str(iter)+".mat","feSolution", "solutionType", "Error", "NumElements", "ElementDegree", "f", "g", "h");
    fprintf(fOut,"%s,%.16E,%d,%d,%s,%s,%s\n",solutionType,Error,NumElements,ElementDegree,string(f),string(g),string(h));
end

fclose(fOut);
