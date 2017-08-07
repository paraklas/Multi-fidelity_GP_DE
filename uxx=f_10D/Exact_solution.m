function u=Exact_solution(x)
    % u = prod(sin(pi*x(1:2:end)),2);
    u = sin(2*pi*x(:,1)).*sin(2*pi*x(:,3));
end