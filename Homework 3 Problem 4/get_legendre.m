function [y_roots,w_values] = get_legendre(n)
    N = n;
    syms y;
    y_roots = vpasolve(legendreP(N,y) == 0);
    temp = diff(legendreP(N,y));
    w_values = 2./((1-y_roots.^2).*subs(temp,y_roots).^2);
    %
    y_roots = double(y_roots);
    w_values = double(w_values);
end