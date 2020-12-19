function F = computeAffineMapping(TargetSpace, Preimage, Domain)
A = [1 Preimage(1); 1 Preimage(2)];
b = [TargetSpace(1); TargetSpace(2)];
c = A\b;

F(Domain) = c(1) + c(2)*Domain;
end
