function f=RHS_H(x)
  %[~, d] = size(x);
  %f = -d*pi^2*Exact_solution(x);
  f = -2*4*pi^2*Exact_solution(x);
end