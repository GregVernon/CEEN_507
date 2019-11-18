%% Exact Solution
clear
close all


syms f(x)

base = sym(0.005);
height = sym(0.005);
I = base*height^3 / 12;
E = sym(1000000);

f(x) = 10*height^3;
g = sym(0);
h = sym(0);
m = sym(0);
q = sym(0);

exactSolution = computeExactSolution(E,I,f,g,h,m,q);
hold on
fplot(exactSolution.U,exactSolution.domain)

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
base = sym(0.005);
height = sym(0.005);
I = base*height^3 / 12;
E = sym(1000000);

f(x) = 10*height^3;
g = sym(0);
h = sym(0);
m = sym(0);
q = sym(0);


nELEM = 2.^(0:7);
for d = 2:3
    eDegree = d;
    for continuity = 1:d-1
        fName = "P" + num2str(d) + "C" + num2str(continuity) + ".csv";
        fOut = fopen(fName,'w+');
        fprintf(fOut,"%s\n","solutionType, Error, y(1),  NumElements, ElementDegree, Continuity, f, g, h, m, q, E, I");
        for n = 1:length(nELEM)
            iter = iter+1;
            disp("Iteration " + num2str(n) + " of " + num2str(length(nELEM)) + "; Number of Elements = " + num2str(nELEM(n)))
            
            contVector = [-1 continuity*ones(1,nELEM(n)-1) -1];
            feSolution = main(nELEM(n), eDegree, contVector, f, g, h, m, q, E*I);
            %         fplot(feSolution.U,[0 1])
            err = computeError(exactSolution, feSolution, "Exact");
            
            solutionType = "Finite Element - Exact Integration";
            Error = err;
            NumElements =  nELEM(n);
            ElementDegree = eDegree;
            save("approx_"+num2str(iter)+".mat","feSolution", "solutionType", "Error", "NumElements", "ElementDegree", "continuity", "f", "g", "h", "m", "q", "E", "I");
            fprintf(fOut,"%s,%.16E,%.16E,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s\n",solutionType,Error,feSolution.U(1),NumElements,ElementDegree,continuity,string(f),string(g),string(h),string(m),string(q),string(E),string(I));
        end
        fclose("all");
    end
end
