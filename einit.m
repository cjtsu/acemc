% compute the corona onset electric field
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function [coronainit, ginit, einit_west, eginit_west] = einit( apt, diameter, ground_diameter, roughness)
% adjust the unit
kt = (diameter < 0.1);
diameter(kt) = diameter(kt) .* 100; 
radius = diameter .* 0.5;
kt = (ground_diameter < 0.1);
ground_diameter(kt) = ground_diameter(kt) .* 100; 
ground_radius = ground_diameter .* 0.5;
delta = (1-0.0065 * apt /293)^4.26;
kt = 30.3 * roughness * sqrt(delta);
coronainit = kt.*( 1 + 0.3 ./ sqrt(radius));
ginit = kt.*( 1 + 0.3 ./sqrt( ground_radius));
einit_west = kt.*(1 + 0.3 ./sqrt(radius * delta));
eginit_west = kt.*(1 + 0.3 ./sqrt(ground_radius * delta));