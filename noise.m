% compute the audible noise. 
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function pdb = noise( gmax,linepos,  ndiv, dsubline, xt, ht)
if ( nargin == 5) ht = 1.; end
m = length( xt);
pdb = zeros(m,1);
if ( ndiv <= 4)
    display('warning: BPA formular cannot be used when the number of subconductors less than 4');
end

deq = 0.58 .* ndiv.^0.48 * dsubline * 1e3; % m ->mm
pwl = -164.6 + 120 .* log10( gmax) + 55 .* log10(deq);
for j = 1:m
    r =  (linepos(:,1) - xt(j)).^2 + ( ht - linepos(:,2)).^2;
    r = sqrt(r);
    sla = 10.^(( pwl - 11.4 .* log10(r) - 5.8)./10);
    pdb(j) = 10 * log10( sum(sla));
end