function u = cantileverDeflection(h, b, L, E)
I = (b * h^3) / 12;
f = 10 * h^3;
u = (f * L^4) / (8 * E * I);
end