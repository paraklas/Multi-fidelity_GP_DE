function [xx] = GetBoundary(dim, n)


for i=1:dim
   bb = sprintf('x%d = linspace(%f, %f, %d);', i, 0.0, 1.0, n(i));
   eval(bb);
end


buf_left = sprintf('[X%d', 1);
buf_right = sprintf('ndgrid(x1', 1);
for i = 2:dim
   buf_left = sprintf('%s, X%d',buf_left,i);
   buf_right = sprintf('%s, x%d',buf_right, i);
end
buf_left = sprintf('%s]',buf_left);
buf_right = sprintf('%s)',buf_right);

buf = sprintf('%s=%s;', buf_left, buf_right);
eval(buf);

xx = zeros(prod(n), dim);
for i = 1:dim
    buf = sprintf('reshape(X%d, prod(n), 1);',i);
    xx(:,i) = eval(buf);
end

idx = (xx<1 & xx >0);

ss = sum(idx,2);

idx = (ss < dim);

xx = xx(idx,:);

% figure(1)
% clf
% if dim == 3
%     plot3(xx(:,1), xx(:,2), xx(:,3),'o');
% elseif dim == 2
%     scatter(xx(:,1), xx(:,2),'o');
% end
% xlim([-1 2]);
% ylim([-1 2]);
% zlim([-1 2]);

end