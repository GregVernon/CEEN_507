function [T, iT] = basisTransform(inputBasis, outputBasis)
T = inputBasis.basis ./ outputBasis.basis;
iT = outputBasis.basis ./ inputBasis.basis;
end