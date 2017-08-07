function f=RHS_H(x)
  [~, d] = size(x);
  f = -d*4*pi^2*Exact_solution(x);
end