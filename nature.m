% compute the natural transmisson power. 
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn


function [imp, pnature] = nature( subconductor, diameter, om, ground_wire, ground_diameter, ground_om, rouearth, total_route, phase_start,sdv)
freq = 50;
[y50,z50]= YZ(subconductor, diameter.*0.5, om, ground_wire,ground_diameter.*0.5,ground_om,freq,rouearth,total_route, phase_start);
total3 = total_route * 3;
YY50 = zeros( total3, total3);
ZZ50 = zeros( total3, total3);
for i = 1:total_route
    t = ((i-1)*3+1):((i-1)*3+3);
    YY50(t,t) = trans(y50(t,t));
    ZZ50(t,t) = trans(z50(t,t));
end

m = exp( sqrt(-1) * 2.0/3.0*pi); % 
a1=[1,1,1;1,m^2,m;1,m,m^2];
A = kron(eye(total_route), a1);
ZZ = A \ ZZ50 * A;  %
YY = A \ YY50 * A;  % 

diagz = diag(imag(ZZ));
diagy = diag(imag(YY));

v = 1:3:(total3-2);
Z11 = diagz(v);
Y11 = diagy(v);
imp = sqrt( Z11./Y11);

pnature = sdv .^2 ./ imp;
pnature = pnature ./ 1e6;    % MVA

return

function T = trans(A)
b = (A(1,1) + A(2,2) + A(3,3))/3;
a = (A(1,2) + A(1,3) + A(2,3))/3;
T = [b a a; a b a; a a b];