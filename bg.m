% compute the b-field at the height of ht. default ht=1
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function b = bg( current, linecentre, sag, xpos, ht)
if ( nargin == 4) ht = 1.; end
yt = ht .* ones( size(xpos));
b = bspace( current, linecentre, sag, xpos, yt);