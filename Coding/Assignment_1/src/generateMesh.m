function [eCONN, x] = generateMesh(xMin,xMax,nElem,degree)
% Create mesh "LM" array per Hughes 1.14
eCONN = zeros(degree+1,nElem);
nCount = 1;
for ii = 1:nElem
    for jj = 1:degree+1
        if jj ~=1
            nCount = nCount + 1;
        end
        eCONN(jj,ii) = nCount;
    end
end

% Create equal-spaced mesh nodes
nNodes = nCount;
x = linspace(xMin,xMax,nNodes);

end
