function f=RHS_L(x)
  f = 0.8*RHS_H(x) - 40*prod(x,2) + 30;
end