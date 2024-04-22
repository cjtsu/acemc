% compute the equallient charge. 
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function [p,q] = qline(u,linecenter,rline)
n = length(u);
x = linecenter(:,1);
y = linecenter(:,2);
d = (2.* y ./ rline).^2;
t1 = repmat(x, 1, n);
t2 = repmat(x.',n, 1);
t1 = t1-t2;    % (x(i)-x(j))
z1 = repmat(y, 1, n);
z2 = repmat(y.', n, 1);
z3 = z1+z2;
z4 = z1-z2;
t1square = t1.^2;  % (x(i)-x(j)).^2
t = ( t1square + z3.^2) ./ ( t1square + z4.^2);
t(1:n+1:end) = d;  %
p = 0.5 .* log(t);
episilon = 1/36/pi*1e-9;
q = p\u.*(2*pi*episilon);

% ---------------- old code ---------------
% m = log(2.* y ./ rline);
% p = diag(m);
% t = zero(n,n);
% for i = 1:n
%     for j = (i+1):n
%          t(i,j) = .5*log(((x(j)-x(i))^2+(y(j)+y(i))^2)/((x(j)-x(i))^2+(y(j)-y(i))^2));
%     end
% end
% p = p + t + t';
%--------------------------------------------