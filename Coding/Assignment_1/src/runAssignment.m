%% Exact Solution
clear
close all

T = table();

x = sym('x','real');
f(x) = x;
g = sym(0);
h = sym(0);

exactSolution = computeExactSolution(f,g,h);
hold on
fplot(exactSolution.U,exactSolution.domain,'k')

T.Solution = {exactSolution};
T.SolutionType = "Exact";
T.Error = 0;
T.NumElements =  inf;
T.f = {f};
T.g = g;
T.h = h;

save("output.mat","T");

%% Test Sweep, # Elements

F = {symfun(1,x), symfun(x,x), symfun(x^2,x)};
for loadCase = 1:length(F)
    f = F{loadCase};
    nELEM = 2.^(1:8);
    for n = 1:length(nELEM)
        disp("Iteration " + num2str(n) + " of " + num2str(length(nELEM)) + "; Number of Elements = " + num2str(nELEM(n)))
        feSolution = main(nELEM(n),2,f,g,h);
        fplot(feSolution.U,[0 1])
        err = computeError(exactSolution, feSolution, "Exact");
        
        TNew = table();
        TNew.Solution = {feSolution};
        TNew.SolutionType = "Finite Element - Exact Integration";
        TNew.Error = err;
        TNew.NumElements =  nELEM(n);
        TNew.f = {f};
        TNew.g = g;
        TNew.h = h;
        T = [T;TNew];
        save("output.mat","T");
    end
end

%% Analysis & Plotting
