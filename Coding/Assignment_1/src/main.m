function feSolution = main(nElem, elemDegree, loadCase)
xMin = 0;
xMax = 1;

[eCONN,nodes] = generateMesh(xMin,xMax,nElem,elemDegree);
bFun = basisFunction("Lagrange", elemDegree, sym("x","real"), [-1 1]);
ELEM = createElements(eCONN,nodes,bFun);

feSolution.Elements = ELEM;
feSolution.Mesh.Connectivity = eCONN;
feSolution.Mesh.Nodes = nodes;
end