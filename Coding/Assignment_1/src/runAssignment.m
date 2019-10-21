%% Exact Solution
clear
close all


x = sym('x','real');
F = {sym(1), x, x^2};
g = sym(0);
h = sym(0);
for loadCase = 1:length(F)
    f = F{loadCase};
    
    exactSolution = computeExactSolution(f,g,h);
    hold on
    fplot(exactSolution.U,exactSolution.domain)
    
    solutionType = "Exact";
    Error = 0;
    NumElements =  inf;
    ElementDegree = NaN;
    
    save("exact_"+num2str(loadCase)+".mat",'exactSolution', 'solutionType', 'Error', 'NumElements', 'ElementDegree', 'f', 'g', 'h');
    exSolution{loadCase} = exactSolution;
end
exactSolution = exSolution;

%% Test Sweep, # Elements
clearvars -except exactSolution
fOut = fopen("TestSweep.csv",'w+');
fprintf(fOut,"%s\n","solutionType, Error, NumElements, ElementDegree, f, g, h");
iter = 0;
x = sym('x','real');
F = {symfun(1,x), symfun(x,x), symfun(x^2,x)};
g = sym(0);
h = sym(0);
for loadCase = 1:length(F)
    f = F{loadCase};
    nELEM = 2.^(1:4);
    eDegree = 3;
    for n = 1:length(nELEM)
        iter = iter+1;
        disp("Iteration " + num2str(n) + " of " + num2str(length(nELEM)) + "; Number of Elements = " + num2str(nELEM(n)))
        feSolution = main(nELEM(n),eDegree,f,g,h);
        %         fplot(feSolution.U,[0 1])
        err = computeError(exactSolution{loadCase}, feSolution, "Exact");
   
        solutionType = "Finite Element - Exact Integration";
        Error = err;
        NumElements =  nELEM(n);
        ElementDegree = eDegree;
        save("approx_"+num2str(iter)+".mat","feSolution", "solutionType", "Error", "NumElements", "ElementDegree", "f", "g", "h");
        fprintf(fOut,"%s,%.16E,%d,%d,%s,%s,%s\n",solutionType,Error,NumElements,ElementDegree,string(f),string(g),string(h));
    end
end
fclose(fOut);

%% Analysis & Plotting
% clear
% Exact = load('exact_2.mat');
% n = 0;
% for ii = 9:16
%     n = n+1;
%     Approx{n,1} = matfile("approx_"+num2str(ii)+".mat");
% end
% 
% % Plot Error vs nElem
% for ii = 1:8
%     err(ii) = Approx{ii}.Error;
%     nElem(ii) = Approx{ii}.NumElements;
% end
% loglog(nElem,err,'k')
% 
% 
% %% Load Case 2
% clear
% Exact = load('exact_2.mat');
% n = 0;
% for ii = 17:24
%     n = n+1;
%     Approx{n,1} = matfile("approx_"+num2str(ii)+".mat");
% end
% 
% % Plot Error vs nElem
% for ii = 1:8
%     err(ii) = Approx{ii}.Error;
%     nElem(ii) = Approx{ii}.NumElements;
% end
% loglog(nElem,err,'k')
