% compute the audible noise using formula in 345 kv handbook
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function pdb = noise345( emax, linepos, ndiv, rsubline, xt, ht)
if ( nargin < 6) ht = 1.; end
n = length( emax); % 
m = length( xt);   % 
rsubline = rsubline .* 1e3;  % m->mm
A1 = 46.4 - 665 ./ emax;
A1 = 10.^(A1./10);
kn = ones(n,1); kn(ndiv==1) = 5.6; kn(ndiv==2) = 1.8;
dsubline = rsubline .* 2;
A = ndiv .^2 .* (dsubline ./38).^4.4 .* A1 .* kn;
A = repmat( A, 1, m);

ndivmod = ndiv;
ndivmod(ndiv<4) = 4;
ec = (1.25 .* dsubline - 4.57)./ ( dsubline ./ 10 - 1.07);
ec = ec ./(1 + 0.027.*( ndivmod - 4));

emec = emax ./ec;
xmc = 10.* ( emec - 0.8);
C = ( 63.4 .* xmc.^2 + 1.87 * xmc.^3 -1.15 * xmc.^4)./1000;
C = repmat(C,1,m);

x1 = repmat( linepos(:,1),1, m);
y1 = repmat( linepos(:,2),1, m);
xt = repmat( xt.', n, 1);
yt = ones(n,m) .* ht;
d = (x1-xt).^2 + (y1-yt).^2;
r = sqrt(d);
j = exp(-0.0075 .* r) .* A ./ 4 ./r;
j = j .* C;
p = 20.5 .* sqrt( sum(j)) ./ 1000;
pdb = 20 .* log10( p ./ 2 .* 1e5);
pdb = pdb';