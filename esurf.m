% compute the surface e-field using the multiple mirror method. 
% accurately enough if we use only 1 more mirror
% part of a package for calculating the electromagtic enviroment of AC power transsion lines.
% initially created around March 2007. 
% optimized using vector around 2009.
% open sourced in 2024.
% code written by chijie@tsinghua.edu.cn

function [ef, ground_max, ...
    conductor_phase_max, conductor_phase_mean]...
        = esurf( q, linecentre, rline, totalroute, phaseinfo)
SQRTJ = sqrt(-1);
episilon2pi = 1/36/pi*1e-9*2*pi;
wireno = length( linecentre);
wirepos = linecentre(:,1) + SQRTJ * linecentre(:,2);
wire1 = repmat( wirepos, 1, wireno);
wire2 = repmat( wirepos.', wireno, 1);  % transpose of wirepos
dif = wire1 - wire2;
wirefi = angle(dif);
wiredis = abs(dif);
q_div_r_mul_2 = 2.* q ./ rline;
q_div_r_mul_2 = q_div_r_mul_2 .'; 
m1 = repmat( q_div_r_mul_2, 360, 1);  
alpha = (0:359) ./360 * 2 * pi; 
tmpalpha = repmat( alpha .', 1, wireno);
lenalpha = length(alpha);
ef = ones(wireno,1);

v = 1:wireno:(wireno*wireno);
wiredis(v) = 1.;
q_div_rline = q ./ rline;

for i = 1:wireno % each subconductor.
    r_div_d = rline(i) ./ wiredis(i,:);
   
    r_div_d(i) = 0;
    r_div_d2 = r_div_d .* r_div_d;
    % r_div_d3 = r_div_d .* r_div_d2;    
    r_div_d = repmat( r_div_d, lenalpha, 1);
    r_div_d2 = repmat( r_div_d2, lenalpha, 1);
    % r_div_d3 = repmat( r_div_d3, lenalpha, 1);     
    tmpwirefi = repmat( wirefi(i,:), lenalpha, 1);
    tmpang = tmpalpha - tmpwirefi;
    costmp = cos(tmpang);
    cos2tmp = 2. * costmp .^2 -1; 
    % cos2tmp = cos(2.*tmpang);
    % cos3tmp = cos(3.*tmpang);      
    % ep = m1.* (r_div_d .* costmp + r_div_d2 .* cos2tmp + r_div_d3 .* cos3tmp);
    ep = m1.* (r_div_d .* costmp + r_div_d2 .* cos2tmp);
    ep = sum(ep,2);
    ep = q_div_rline(i) - ep;    
    ep = abs(ep);
    ef(i) = max(ep);   
end
ef = ef ./ episilon2pi;
cnt = totalroute * 3;  % 
groundcnt = wireno - phaseinfo(totalroute+1, 1) + 1;  % 
ground_max = ef(1:groundcnt);
ph = reshape( phaseinfo.', [], 1);
ph = ph + groundcnt;   % linecenter
conductor_phase_max = zeros(cnt, 1);
conductor_phase_mean = zeros(cnt,1);
for i=1:cnt
    tf = ef(ph(i):ph(i+1)-1);
    conductor_phase_max(i) = max(tf);
    conductor_phase_mean(i) = mean(tf);
end
ef = ef ./1e5;
ground_max  = ground_max ./ 1e5;
conductor_phase_max = conductor_phase_max ./ 1e5;
conductor_phase_mean = conductor_phase_mean ./ 1e5;