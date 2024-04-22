% compute the e-field at (xposition, yposition)
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function E = espace( q, linepos, sag, xposition, yposition)
n = length(xposition);
m = length(q);
x = linepos(:,1);
y = linepos(:,2) - sag ./ 3;  % consider the sag
x1 = repmat( x, 1, n);
y1 = repmat( y, 1, n);
qn = repmat( q, 1, n);
% assume xpos is a row vector
if ( size(xposition, 1) == 1) 
    x2 = repmat( xposition, m, 1);
else
    x2 = repmat( xposition', m, 1);
end
if ( size(yposition, 1) == 1) 
    y2 = repmat( yposition, m, 1);
else
    y2 = repmat( yposition', m, 1);
end
kx2 = (x1-x2).^2;
l1 = kx2 + ( y1-y2).^2;
l2 = kx2 + ( y1+y2).^2;
ex = qn .* ( (x2-x1) ./ l1 - (x2-x1) ./ l2);
ey = qn .* ( (y2-y1) ./ l1 - (y1+y2) ./ l2);
Ex = sum(ex);
Ey = sum(ey);
episilon = 1 / 36 / pi * 1e-9;
Ex = Ex./(2 * pi * episilon);
Ey = Ey./(2 * pi * episilon);
% Ex = Exsin(wt+faix) Ey = Eysin(wt+faiy)
faix = atan(imag(Ex)./real(Ex)); 
faiy = atan(imag(Ey)./real(Ey));
Ex = abs(Ex); Ey = abs(Ey);
% calcualte the maximum of abs(Ex^2+Ey^2) analyticlly
a = Ex.^2; b = Ey.^2;
c = a .* cos(2.*faix) + b .* cos(2.*faiy);
d = a .* sin(2.*faix) + b .* sin(2.*faiy);
p = ( a+b + sqrt(c.^2+d.^2))./2; 
E = sqrt(p);