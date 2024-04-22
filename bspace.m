% compute the b-field at (xpos, ypos)
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function bs = bspace( current, linepos, sag, xpos, ypos)
n = length( xpos);
m = length( current);
mu = 2*10^(-7);
x = linepos(:,1);
 % assume the line is with sag, and the averged y should be lower.
 % averge height of the line = y-sag/3
y = linepos(:,2) - sag ./ 3; 
x1 = repmat( x, 1, n);
y1 = repmat( y, 1, n);
cn = repmat( current, 1, n);
% assume xpos is a row vector
if ( size(xpos, 1) == 1) 
    x2 = repmat( xpos, m, 1);
else
    x2 = repmat( xpos', m, 1);
end
if ( size(ypos, 1) == 1) 
    y2 = repmat( ypos, m, 1);
else
    y2 = repmat( ypos', m, 1);
end
l1 = (x1-x2).^2 + (y1-y2).^2;
bx = cn .* ( x1-x2) ./l1;
by = cn .* ( y1-y2) ./l1;
bx = sum(bx); by = sum(by);
faix = angle(bx);
faiy = angle(by);
bx = abs(bx); by = abs(by);
% B = i*Bx+j*By, and should output max(abs(B))
% caluclte the max(abs(B)) by analytic solution
a = bx.^2; b = by.^2;
c = a .* cos(2.*faix) + b .* cos(2.*faiy);
d = a .* sin(2.*faix) + b .* sin(2.*faiy);
p = (a+b+ sqrt(c.^2+d.^2))./2; 
bs = sqrt(p) .* mu;  