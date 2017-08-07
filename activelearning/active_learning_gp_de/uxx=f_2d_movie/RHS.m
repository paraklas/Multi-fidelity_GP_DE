function f=RHS(x)
  [~, d] = size(x);
  %f = -d*4*pi^2*Exact_solution(x);
  f = -d*pi^2*Exact_solution(x);
end