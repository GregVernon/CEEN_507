function u = cantileverDeflection(h, b, L, E)
I = (b * h^3) / 12;
f = 10 * h^3;
u = (f * L^4) / (8 * E * I);
end

function u(x) = cantileverXDeflection(h, b, L, E)
I = (b * h^3) / 12;
f = 10 * h^3;
u(x) = (f * x^2) / (24 * E * I) * (x^2 + 6 * L^2 - 4*L*x);
end
