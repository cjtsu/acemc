% compute the e-field at ht
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function E = eg( q, linepos, sag, xt, ht)
if ( nargin == 4) ht = 1.; end
yt = ht .* ones( size(xt));
E = espace( q, linepos, sag, xt, yt);