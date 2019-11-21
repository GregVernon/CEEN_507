%% Exact Solution
clear
close all


syms f(x) g h
f(x) = x^2;
g = sym(0);
h = sym(0);

exactSolution = computeExactSolution(f,g,h);
% hold on
% fplot(exactSolution.U,exactSolution.domain)

solutionType = "Exact";
Error = 0;
NumElements =  inf;
ElementDegree = NaN;

save("exact.mat",'exactSolution', 'solutionType', 'Error', 'NumElements', 'ElementDegree', 'f', 'g', 'h');
exSolution = exactSolution;
exactSolution = exSolution;

%% Test Sweep, # Elements
clearvars -except exactSolution
iter = 0;

x = sym('x','real');
f(x) = x^2;
g = sym(0);
h = sym(0);


nELEM = 2.^(0:7);
for p = 2:3
    eDegree = p;
    for continuity = 1:p-1
        fName = "P" + num2str(p) + "C" + num2str(continuity) + ".csv";
        fOut = fopen(fName,'w+');
        fprintf(fOut,"%s\n","solutionType, Error u(x), Error du(x), Error d2u(x), Error d3u(x), y(1),  NumElements, ElementDegree, Continuity, ndofs, nnodes, f, g, h");
        for n = 1:length(nELEM)
            iter = iter+1;
            disp("Iteration " + num2str(n) + " of " + num2str(length(nELEM)) + "; Number of Elements = " + num2str(nELEM(n)))
            
            elemContinuity = [-1 continuity*ones(1,nELEM(n)-1) -1];
            feSolution = main(nELEM(n),eDegree,elemContinuity,f,g,h);
            %         fplot(feSolution.U,[0 1])
            err = zeros(4,1);
            for deriv = 0:3
                err(deriv+1) = computeError(exactSolution, feSolution, "Exact",deriv);
            end
            
            solutionType = "Finite Element - Exact Integration";
            Error = double(err);
            NumElements =  nELEM(n);
            ElementDegree = eDegree;
			ndofs = length(feSolution.LinearSystem.F);
			nnodes = length(feSolution.LinearSystem.d);
			
            save("approx_"+num2str(iter)+".mat","feSolution", "solutionType", "Error", "NumElements", "ElementDegree","continuity", "f", "g", "h");
            fprintf(fOut,"%s,%.16E,%.16E,%.16E,%.16E,%.16E,%d,%d,%d,%d,%d,%s,%s,%s\n",solutionType,Error(1),Error(2),Error(3),Error(4),feSolution.U(0),NumElements,ElementDegree,continuity,ndofs,nnodes,string(f),string(g),string(h));
        end
        fclose('all');
    end
end
