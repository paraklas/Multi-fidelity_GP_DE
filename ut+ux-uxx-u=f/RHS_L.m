function f=RHS_L(x)
    f = 0.8*RHS_H(x) - 5*prod(x,2) - 20;
end