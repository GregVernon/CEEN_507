clear

splineSpace.degree = 4;
splineSpace.nodes = [0 1 2 3 4 5];
splineSpace.continuityVector = [-1 3 2 1 0 -1];

b = bspline(splineSpace);
% b.basis.functions

% fplot(b.basis.functions,[-1 6])
% axis equal

[b, T] = b.bezierExtraction();
