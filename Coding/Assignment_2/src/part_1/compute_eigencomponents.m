function I = compute_eigencomponents(x, eigvec)
I = zeros(length(eigvec),1);
eigvec = real(double(eigvec));
x = real(double(x));
for ii = 1:size(eigvec,2)
    I(ii) = dot(eigvec(:,ii),x);
end
end