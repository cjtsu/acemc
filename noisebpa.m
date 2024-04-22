% compute the audible noise using bpa formula. 
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function pdb = noisebpa( emax,linepos, ndiv, rsubline, xt, ht)
rsubline = rsubline .* 1e3;   % m ->mm
if ( nargin == 5) ht = 1.; end
m = length( xt);
n = length(linepos(:,1));
if ( min(ndiv) <= 4)
%    display('warning: BPA formular cannot be used when the number of subconductors less than 4');
end
dsubline = rsubline .* 2;
deq = 0.58 .* ndiv.^0.48 .* dsubline;  %
pwl = -164.6 + 120 .* log10(emax) + 55 .* log10(deq);
reppwl = repmat( pwl, 1, m);
x1 = repmat( linepos(:,1), 1, m);
y1 = repmat( linepos(:,2), 1, m);
y2 = ht .* ones(n, m);
x2 = repmat( xt.', n, 1);
d = (x1-x2).^2 + (y1-y2).^2;
sla = 10.^((reppwl -11.4 * 0.5.*log10(d) -5.8)./10);  
pdb = 10 .*log10(sum(sla,1));
pdb = pdb';