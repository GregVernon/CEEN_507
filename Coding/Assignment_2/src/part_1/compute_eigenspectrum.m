function [eigvec,eigval] = compute_eigenspectrum(A)
[eigvec,eigval] = eig(full(A));
eigval = diag(eigval);
for ii = 1:length(eigval)
    eigvec(:,ii) = eigvec(:,ii) ./ norm(eigvec(:,ii));
end

[~,sidx] = sort(real(double(eigval)),"ascend");
eigval = eigval(sidx);
eigvec = eigvec(:,sidx);
end