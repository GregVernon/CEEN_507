%% Mesh Generation
% Create mesh "LM" array per Hughes 1.14
[eCONN,x] = generateMesh(0,2,2,1);
assert(isequal(eCONN,[1 2; 2 3]));
assert(isequal(x,[0 1 2]));

[eCONN,x] = generateMesh(0,2,2,2);
assert(isequal(eCONN,[1 2 3; 3 4 5]'));
assert(isequal(x,[0 0.5 1 1.5 2]));

[eCONN,x] = generateMesh(0,3,3,1);
assert(isequal(eCONN,[1 2; 2 3; 3 4 ]'));
assert(isequal(x,[0 1 2 3]));

[eCONN,x] = generateMesh(0,3,3,2);
assert(isequal(eCONN,[1 2 3; 3 4 5; 5 6 7]'));
assert(isequal(x,[0 0.5 1 1.5 2 2.5 3]));
